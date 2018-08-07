#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
	Kover: Learn interpretable computational phenotyping models from k-merized genomic data
	Copyright (C) 2018  Alexandre Drouin & Gael Letarte

	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import gc
import h5py as h
import logging
import numpy as np
import pandas as pd

from math import ceil
from os import getpid, mkdir, listdir
from os.path import basename, exists, getsize, join, splitext
from shutil import rmtree
from time import time
from uuid import uuid1

from ..utils import _minimum_uint_size, _pack_binary_bytes_to_ints
from .tools.kmer_count import contigs_count_kmers, reads_count_kmers
from .tools.kmer_pack import contigs_pack_kmers, reads_pack_kmers

KMER_MATRIX_PACKING_SIZE = 64
KMER_MATRIX_DTYPE = np.uint64
PHENOTYPE_LABEL_DTYPE = np.uint8
BLOCK_SIZE = 100000


def _create_hdf5_file_no_chunk_caching(path):
    # Create the HDF5 File
    h5py_file = h.File(path, "w")

    # Get the default cache properties
    access_property_list = h5py_file.fid.get_access_plist()
    cache_properties = list(access_property_list.get_cache())

    h5py_file.close()

    # Disable chunk caching
    cache_properties[2] = 0  # No chunk caching
    access_property_list.set_cache(*cache_properties)
    file_id = h.h5f.open(path, h.h5f.ACC_RDWR, fapl=access_property_list)

    # Reopen the file without a chunk cache
    h5py_file = h.File(file_id)

    return h5py_file


def _parse_metadata(metadata_path, matrix_genome_ids, warning_callback, error_callback):
    """
	Parses metadata (genome_id{tab}label)
	"""
    logging.debug("Parsing metadata.")
    md_genome_ids, md_genome_labels = zip(*(l.split() for l in open(metadata_path, "r")))
    md_unique_labels, indices = np.unique(md_genome_labels, return_inverse=True)

    # Check case label are 0 and 1 (To be backward compatible with Kover previous dataset creation version)
    if not(len(md_unique_labels) == 2 and '0' in md_unique_labels and '1' in md_unique_labels):
        # Sorting labels alphabetically for consistent indices assignement across multiple datasets
        md_unique_labels.sort()
        label_to_indice = {md_unique_labels[l]:l for l in range(len(md_unique_labels))}
        indices = np.array([label_to_indice[l] for l in md_genome_labels])

    if len(md_unique_labels) < 2:
        error_callback(Exception("The dataset must contain at least 2 different phenotypes"))

    elif len(md_unique_labels) > 255:
        error_callback(Exception("The dataset can contain at most 255 different phenotypes"))

    elif len(md_unique_labels) == 2:
        classification_type = "binary"

    else:
        classification_type = "multiclass"
    logging.debug("The dataset problem type is " + classification_type + " classification.")

    # Converting labels to numerical ascending values
    numerical_labels = np.arange(0, len(md_unique_labels))
    md_genome_labels = numerical_labels[indices]

    if len(md_genome_ids) > len(set(md_genome_ids)):
        error_callback(Exception("The metadata contains multiple values for the same genome."))

    genomes_only_in_matrix = set(matrix_genome_ids) - set(md_genome_ids)
    if len(genomes_only_in_matrix) > 0:
        warning_callback("Missing metadata for %d genomes (%s). These genomes will be discarded." % (
            len(genomes_only_in_matrix), ", ".join(genomes_only_in_matrix)))
    del genomes_only_in_matrix

    genomes_only_in_metadata = set(md_genome_ids) - set(matrix_genome_ids)
    if len(genomes_only_in_metadata) > 0:
        warning_callback("The metadata contains values for %d genomes that are not in the genomic data (%s)." % (
            len(genomes_only_in_metadata), ", ".join(genomes_only_in_metadata)))
    del genomes_only_in_metadata

    keep_genome_ids, keep_genome_labels = zip(
        *((md_genome_ids[i], md_genome_labels[i]) for i in xrange(len(md_genome_ids)) if
          md_genome_ids[i] in matrix_genome_ids))

    return np.array(keep_genome_ids), np.array(keep_genome_labels, dtype=np.uint8), md_unique_labels, classification_type


def from_tsv(tsv_path, output_path, phenotype_description, phenotype_metadata_path, gzip, warning_callback=None,
             error_callback=None, progress_callback=None):
    def get_kmer_length(tsv_path):
        with open(tsv_path, "r") as f:
            f.next()
            kmer_len = len(f.next().split("\t")[0])
        return kmer_len

    def get_kmer_count(tsv):
        # XXX: This is can break if the file format changes!
        #      Assumes that each line has the format kmer_seq\tV\tV\t...V\n with one binary V per genome
        total_size = getsize(tsv)
        with open(tsv, "r") as f:
            header_size = len(f.next())
            line_size = len(f.next())
        content_size = float(total_size) - header_size
        if content_size % line_size != 0:
            raise Exception()
        return int(content_size / line_size)

    # Execution callback functions
    if warning_callback is None:
        warning_callback = lambda w: logging.warning(w)
    if error_callback is None:
        def normal_raise(exception):
            raise exception

        error_callback = normal_raise
    if progress_callback is None:
        progress_callback = lambda t, p: None

    if (phenotype_description is None and phenotype_metadata_path is not None) or (
                    phenotype_description is not None and phenotype_metadata_path is None):
        raise ValueError("If a phenotype is specified, it must have a description and a metadata file.")

    kmer_len = get_kmer_length(tsv_path)
    kmer_count = get_kmer_count(tsv_path)

    kmer_dtype = 'S' + str(kmer_len)
    kmer_by_matrix_column_dtype = _minimum_uint_size(kmer_count)
    compression = "gzip" if gzip > 0 else None
    compression_opts = gzip if gzip > 0 else None
    tsv_block_size = min(kmer_count, BLOCK_SIZE)

    h5py_file = _create_hdf5_file_no_chunk_caching(output_path)
    h5py_file.attrs["created"] = time()
    h5py_file.attrs["uuid"] = str(uuid1())
    h5py_file.attrs["genome_source_type"] = "tsv"
    h5py_file.attrs["genomic_data"] = tsv_path
    h5py_file.attrs["phenotype_description"] = phenotype_description if phenotype_description is not None else "NA"
    h5py_file.attrs[
        "phenotype_metadata_source"] = phenotype_metadata_path if phenotype_metadata_path is not None else "NA"
    h5py_file.attrs["compression"] = "gzip (level %d)" % gzip

    # Read list of genome identifiers
    reader = pd.read_table(tsv_path, sep='\t', index_col=0, iterator=True, engine="c")
    genome_ids = reader.get_chunk(1).columns.values
    del reader
    logging.debug("The k-mer matrix contains %d genomes." % len(genome_ids))
    if len(set(genome_ids)) < len(genome_ids):
        error_callback(Exception("The genomic data contains genomes with the same identifier."))

    # Extract/write the metadata
    if phenotype_description is not None:
        genome_ids, labels,\
        labels_tags, classification_type = _parse_metadata(metadata_path=phenotype_metadata_path,
                                                           matrix_genome_ids=genome_ids,
                                                           warning_callback=warning_callback,
                                                           error_callback=error_callback)

        h5py_file.attrs["classification_type"] = classification_type
        # Sort the genomes by label for optimal better performance
        logging.debug("Sorting genomes by metadata label for optimal performance.")
        sorter = np.argsort(labels)
        genome_ids = genome_ids[sorter]
        labels = labels[sorter]
        logging.debug("Creating the phenotype metadata dataset.")
        phenotype = h5py_file.create_dataset("phenotype", data=labels, dtype=PHENOTYPE_LABEL_DTYPE)
        phenotype.attrs["description"] = phenotype_description
        del phenotype, labels

    # Write genome ids
    logging.debug("Creating the genome identifier dataset.")
    h5py_file.create_dataset("genome_identifiers",
                             data=genome_ids,
                             compression=compression,
                             compression_opts=compression_opts)

    # Write labels tags
    logging.debug("Creating the phenotype tags dataset.")
    h5py_file.create_dataset("phenotype_tags",
                             data=labels_tags,
                             compression=compression,
                             compression_opts=compression_opts)

    # Initialize kmers (kmer_list) dataset
    logging.debug("Creating the kmer sequence dataset.")
    kmers = h5py_file.create_dataset("kmer_sequences",
                                     shape=(kmer_count,),
                                     dtype=kmer_dtype,
                                     compression=compression,
                                     compression_opts=compression_opts)

    # Initialize kmer_matrix dataset
    logging.debug("Creating the kmer matrix dataset.")
    kmer_matrix = h5py_file.create_dataset("kmer_matrix",
                                           shape=(
                                               int(ceil(1.0 * len(genome_ids) / KMER_MATRIX_PACKING_SIZE)), kmer_count),
                                           dtype=KMER_MATRIX_DTYPE,
                                           compression=compression,
                                           compression_opts=compression_opts,
                                           chunks=(1, tsv_block_size))

    # Initialize kmer_by_matrix_column dataset
    logging.debug("Creating the kmer sequence/matrix column mapping dataset.")
    kmer_by_matrix_column = h5py_file.create_dataset("kmer_by_matrix_column",
                                                     shape=(kmer_count,),
                                                     dtype=kmer_by_matrix_column_dtype,
                                                     compression=compression,
                                                     compression_opts=compression_opts)

    logging.debug("Transferring the data from TSV to HDF5.")
    tsv_reader = pd.read_table(tsv_path, index_col='kmers', sep='\t', chunksize=tsv_block_size)
    n_blocks = int(ceil(1.0 * kmer_count / tsv_block_size))
    n_copied_blocks = 0.
    for i, chunk in enumerate(tsv_reader):
        logging.debug("Block %d/%d." % (i + 1, n_blocks))
        progress_callback("Creating", n_copied_blocks / n_blocks)
        logging.debug("Reading data from TSV file.")
        kmers_data = chunk.index.values.astype(kmer_dtype)
        read_block_size = kmers_data.shape[0]
        block_start = i * tsv_block_size
        block_stop = block_start + read_block_size
        n_copied_blocks += 0.5
        progress_callback("Creating", n_copied_blocks / n_blocks)

        if block_start > kmer_count:
            break

        logging.debug("Writing data to HDF5.")
        kmers[block_start:block_stop] = kmers_data
        attribute_classification_sorted_by_strains = chunk[genome_ids]
        attribute_classification_data = attribute_classification_sorted_by_strains.T.values.astype(np.uint8)
        logging.debug("Packing the data.")
        kmer_matrix[:, block_start:block_stop] = _pack_binary_bytes_to_ints(attribute_classification_data,
                                                                            pack_size=KMER_MATRIX_PACKING_SIZE)
        n_copied_blocks += 0.5
        progress_callback("Creating", n_copied_blocks / n_blocks)

        # Write kmer_by_matrix_column
        kmer_by_matrix_column[block_start:block_stop] = np.arange(block_start, block_stop, dtype=kmer_by_matrix_column)
        logging.debug("Garbage collection.")
        gc.collect()  # Clear the memory objects created during the iteration, or else the memory will keep growing.

    h5py_file.close()

    logging.debug("Dataset creation completed.")


def from_contigs(contig_list_path, output_path, kmer_size, filter_singleton, phenotype_description, phenotype_metadata_path,
                 gzip, temp_dir, nb_cores, verbose, progress, warning_callback=None,
                 error_callback=None):
    compression = "gzip" if gzip > 0 else None
    compression_opts = gzip if gzip > 0 else None

    # Execution callback functions
    if warning_callback is None:
        warning_callback = lambda w: logging.warning(w)
    if error_callback is None:
        def normal_raise(exception):
            raise exception
        error_callback = normal_raise

    # Make sure that the tmp data is unique to the current process
    temp_dir = join(temp_dir, str(getpid()))
    if not exists(temp_dir):
        mkdir(temp_dir)

    if (phenotype_description is None and phenotype_metadata_path is not None) or (
                    phenotype_description is not None and phenotype_metadata_path is None):
        error_callback(ValueError("If a phenotype is specified, it must have a description and a metadata file."))

    # Find the contig file for each genome and verify that it exists
    contig_file_by_genome_id = dict(l.split() for l in open(contig_list_path, "r"))
    for g_id, contig_file in contig_file_by_genome_id.iteritems():
        if not exists(contig_file):
            error_callback(IOError("The contig file for genome %s cannot be found: %s" % (str(g_id), contig_file)))

    logging.debug("The k-mer matrix contains %d genomes." % len(contig_file_by_genome_id))
    if len(set(contig_file_by_genome_id.keys())) < len(contig_file_by_genome_id.keys()):
        error_callback(Exception("The genomic data contains genomes with the same identifier."))

    h5py_file = _create_hdf5_file_no_chunk_caching(output_path)
    h5py_file.attrs["created"] = time()
    h5py_file.attrs["uuid"] = str(uuid1())
    h5py_file.attrs["genome_source_type"] = "contigs"
    h5py_file.attrs["genomic_data"] = contig_list_path
    h5py_file.attrs["phenotype_description"] = phenotype_description if phenotype_description is not None else "NA"
    h5py_file.attrs[
        "phenotype_metadata_source"] = phenotype_metadata_path if phenotype_metadata_path is not None else "NA"
    h5py_file.attrs["filter"] = filter_singleton
    h5py_file.attrs["compression"] = "gzip (level %d)" % gzip

    # Extract/write the metadata
    if phenotype_description is not None:
        genome_ids, labels,\
        labels_tags, classification_type = _parse_metadata(metadata_path=phenotype_metadata_path,
                                                           matrix_genome_ids=contig_file_by_genome_id.keys(),
                                                           warning_callback=warning_callback,
                                                           error_callback=error_callback)

        h5py_file.attrs["classification_type"] = classification_type

        # Sort the genomes by label for optimal better performance
        logging.debug("Sorting genomes by metadata label for optimal performance.")
        sorter = np.argsort(labels)
        genome_ids = genome_ids[sorter]
        labels = labels[sorter]
        logging.debug("Creating the phenotype metadata dataset.")
        phenotype = h5py_file.create_dataset("phenotype", data=labels, dtype=PHENOTYPE_LABEL_DTYPE)
        phenotype.attrs["description"] = phenotype_description
        del phenotype, labels

    # Write genome ids
    logging.debug("Creating the genome identifier dataset.")
    h5py_file.create_dataset("genome_identifiers",
                             data=genome_ids,
                             compression=compression,
                             compression_opts=compression_opts)

    # Write labels tags
    logging.debug("Creating the phenotype tags dataset.")
    h5py_file.create_dataset("phenotype_tags",
                             data=labels_tags,
                             compression=compression,
                             compression_opts=compression_opts)

    h5py_file.close()

    logging.debug("Initializing DSK.")

    # Preparing input file for multidsk
    files_sorted = ["%s\n" % contig_file_by_genome_id[id] for id in genome_ids]
    open(join(temp_dir, "list_contigs_files"), "w").writelines(files_sorted)

    # Calling multidsk
    contigs_count_kmers(file_path=join(temp_dir, "list_contigs_files"),
                        out_dir=temp_dir,
                        kmer_size=kmer_size,
                        out_compress=gzip,
                        nb_cores=nb_cores,
                        verbose=int(verbose),
                        progress=progress)
    logging.debug("K-mers counting completed.")

    # Preparing input file for dsk2kover
    list_contigs = [join(temp_dir, basename(splitext(file)[0]) + ".h5") for file in files_sorted]
    file_dsk_output = open(join(temp_dir, "list_h5"), "w")
    for line in list_contigs:
        file_dsk_output.write(line + "\n")
    file_dsk_output.close()

    # Calling dsk2kover
    logging.debug("Initializing DSK2Kover.")
    contigs_pack_kmers(file_path=join(temp_dir, "list_h5"),
                       out_path=output_path,
                       filter_singleton=filter_singleton,
                       kmer_length=kmer_size,
                       compression=gzip,
                       chunk_size=BLOCK_SIZE,
                       nb_genomes=len(genome_ids),
                       progress=progress)

    logging.debug("Removing temporary files.")
    rmtree(temp_dir)

    # progress_callback("dsk2kover", 1)
    logging.debug("Dataset creation completed.")


def from_reads(reads_folders_list_path, output_path, kmer_size, abundance_min, filter_singleton, phenotype_description,
               phenotype_metadata_path, gzip, temp_dir, nb_cores, verbose, progress, warning_callback=None,
                 error_callback=None):
    supported_extensions = ['.fastq','.fastq.gz']
    compression = "gzip" if gzip > 0 else None
    compression_opts = gzip if gzip > 0 else None

    # Execution callback functions
    if warning_callback is None:
        warning_callback = lambda w: logging.warning(w)
    if error_callback is None:
        def normal_raise(exception):
            raise exception
        error_callback = normal_raise

    # Make sure that the tmp data is unique to the current process
    temp_dir = join(temp_dir, str(getpid()))
    if not exists(temp_dir):
        mkdir(temp_dir)

    if (phenotype_description is None and phenotype_metadata_path is not None) or (
                    phenotype_description is not None and phenotype_metadata_path is None):
        error_callback(ValueError("If a phenotype is specified, it must have a description and a metadata file."))

    # Find the read folder for each genome and verify that it exists
    reads_folder_by_genome_id = dict(l.split() for l in open(reads_folders_list_path, "r"))
    for g_id, read_dir in reads_folder_by_genome_id.iteritems():
        if not exists(read_dir):
            error_callback(IOError("The read directory for genome %s cannot be found: %s" % (str(g_id), read_dir)))

    logging.debug("The k-mer matrix contains %d genomes." % len(reads_folder_by_genome_id))
    if len(set(reads_folder_by_genome_id.keys())) < len(reads_folder_by_genome_id.keys()):
        error_callback(Exception("The genomic data contains genomes with the same identifier."))

    h5py_file = _create_hdf5_file_no_chunk_caching(output_path)
    h5py_file.attrs["created"] = time()
    h5py_file.attrs["uuid"] = str(uuid1())
    h5py_file.attrs["genome_source_type"] = "reads"
    h5py_file.attrs["genomic_data"] = reads_folders_list_path
    h5py_file.attrs["phenotype_description"] = phenotype_description if phenotype_description is not None else "NA"
    h5py_file.attrs[
        "phenotype_metadata_source"] = phenotype_metadata_path if phenotype_metadata_path is not None else "NA"
    h5py_file.attrs["filter"] = filter_singleton
    h5py_file.attrs["compression"] = "gzip (level %d)" % gzip

    # Extract/write the metadata
    if phenotype_description is not None:
        genome_ids, labels,\
        labels_tags, classification_type = _parse_metadata(metadata_path=phenotype_metadata_path,
                                                           matrix_genome_ids=reads_folder_by_genome_id.keys(),
                                                           warning_callback=warning_callback,
                                                           error_callback=error_callback)

        h5py_file.attrs["classification_type"] = classification_type
        # Sort the genomes by label for optimal better performance
        logging.debug("Sorting genomes by metadata label for optimal performance.")
        sorter = np.argsort(labels)
        genome_ids = genome_ids[sorter]
        labels = labels[sorter]
        logging.debug("Creating the phenotype metadata dataset.")
        phenotype = h5py_file.create_dataset("phenotype", data=labels, dtype=PHENOTYPE_LABEL_DTYPE)
        phenotype.attrs["description"] = phenotype_description
        del phenotype, labels

    # Write genome ids
    logging.debug("Creating the genome identifier dataset.")
    h5py_file.create_dataset("genome_identifiers",
                             data=genome_ids,
                             compression=compression,
                             compression_opts=compression_opts)

    # Write labels tags
    logging.debug("Creating the phenotype tags dataset.")
    h5py_file.create_dataset("phenotype_tags",
                             data=labels_tags,
                             compression=compression,
                             compression_opts=compression_opts)

    h5py_file.close()

    logging.debug("Initializing DSK.")

    # Preparing input file for multidsk
    files_sorted = []
    list_reads_dsk_output = []
    for id in genome_ids:
        files = [join(reads_folder_by_genome_id[id], file) for file in listdir(reads_folder_by_genome_id[id])
                 if file.endswith(tuple(supported_extensions))]  # Supported extensions
        files_sorted.append(",".join(files) + "\n")
        list_reads_dsk_output.append(join(temp_dir, basename(splitext(files[-1])[0]) + ".h5"))
    open(join(temp_dir, "list_reads_files"), "w").writelines(files_sorted)

    # Calling multidsk
    reads_count_kmers(file_path=join(temp_dir, "list_reads_files"),
                      out_dir=temp_dir,
                      kmer_size=kmer_size,
                      abundance_min=abundance_min,
                      out_compress=gzip,
                      nb_cores=nb_cores,
                      verbose=int(verbose),
                      progress=progress)
    logging.debug("K-mers counting completed.")

    # Preparing input file for dsk2kover
    file_dsk_output = open(join(temp_dir, "list_h5"), "w")
    for line in list_reads_dsk_output:
        file_dsk_output.write(line + "\n")
    file_dsk_output.close()

    # Calling dsk2kover
    logging.debug("Initializing DSK2Kover.")
    reads_pack_kmers(file_path=join(temp_dir, "list_h5"),
                     out_path=output_path,
                     filter_singleton=filter_singleton,
                     kmer_length=kmer_size,
                     compression=gzip,
                     chunk_size=BLOCK_SIZE,
                     nb_genomes=len(genome_ids),
                     progress=progress)

    logging.debug("Removing temporary files.")
    rmtree(temp_dir)

    # progress_callback("dsk2kover", 1)
    logging.debug("Dataset creation completed.")
