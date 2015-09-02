#!/usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import logging
import numpy as np
import pandas as pd
from math import ceil

def _pack_binary_bytes_to_ints(a, pack_size):
    """
    Packs binary values stored in bytes into ints
    """
    logging.debug("Packing block.")
    if pack_size == 64:
        type = np.uint64
    elif pack_size == 32:
        type = np.uint32
    else:
        raise ValueError("Supported data types are 32-bit and 64-bit integers.")

    b = np.zeros((int(ceil(1.0 * a.shape[0] / pack_size)), a.shape[1]), dtype=type)

    packed_rows = 0
    packing_row = 0
    for i in xrange(a.shape[0]):
        if packed_rows == pack_size:
            packed_rows = 0
            packing_row += 1
        tmp = np.asarray(a[i], dtype=type)
        tmp = np.left_shift(tmp, type(pack_size - packed_rows - 1))
        np.bitwise_or(b[packing_row], tmp, out=b[packing_row])
        packed_rows += 1

    logging.debug("Packing done.")

    return b

def get_mic_dataset(dataset_file, strains_to_load):
    f = open(dataset_file, "r")
    strains = []
    sir_labels = []
    for l in f:
        l = l.split()
        strain = l[0]

        # If a strain is not in the list of strains to load, skip it.
        if not strain in strains_to_load:
            logging.debug("Sample %s is in the metadata but not in the k-mer matrix." % strain)
            continue

        strains.append(l[0])
        sir_labels.append(l[1])
    f.close()

    if not np.all(np.unique(strains) != np.array([0, 1])):
        raise RuntimeError("The sample labels must be binary (0, 1) and there must be samples of each class.")

    return np.array(strains), np.array(sir_labels, dtype=np.uint8)

def get_kmer_count_and_len(file):
    f = open(file, "r")

    # Skip the file header
    f.next()

    kmer_len = None
    kmer_count = 0
    for l in f:
        kmer = l.split('\t')[0]
        if kmer_len is None:
            kmer_len = len(kmer)
        else:
            if len(kmer) != kmer_len:
                raise RuntimeError("All k-mers must be of equal length. The expected length was %d and got %d at line %d" % (kmer_len, len(kmer), kmer_count))

        kmer_count += 1

        if kmer_count % 1000 == 0:
            logging.debug("Line %d" % (kmer_count))

    logging.debug("Line %d" % (kmer_count))

    return kmer_count, kmer_len

def convert(kmer_matrix, metadata, output, kmer_len=None, kmer_count=None, gzip=None):
    # Determine the number of k-mers and their length
    if kmer_count is None or kmer_len is None:
        logging.debug("The kmer count or kmer length is unknown. Determining these values by going through the file. (Could be long)")
        kmer_count, kmer_len = get_kmer_count_and_len(kmer_matrix)

    kmer_dtype = 'S' + str(kmer_len)
    attribute_classifications_dtype = np.uint64
    attribute_classifications_packing_size = 64  #TODO: This could be a user specified parameter
    type_by_attribute_dtype = np.uint8
    kmer_by_attribute_dtype = np.uint32  #TODO: This has to be larger than the number of k-mers
    compression = "gzip" if gzip > 0 else None #TODO: The choice of compression method could be user specified
    compression_opts = gzip if gzip > 0 else None
    tsv_block_size = min(kmer_count, 100000)  #TODO: This needs to be smaller or equal to the line length
    #TODO: the tsv block size is currently the same as the chunk size. Is this an optimal choice?
    # TODO: the block size should be a parameter to control memory usage. The default should be fairly small.

    # Create the HDF5 File
    h5py_file = h5py.File(output, "w")

    # Get the default cache properties
    access_property_list = h5py_file.fid.get_access_plist()
    cache_properties = list(access_property_list.get_cache())

    h5py_file.close()

    # Disable chunk caching
    cache_properties[2] = 0  # No chunk caching
    access_property_list.set_cache(*cache_properties)
    file_id = h5py.h5f.open(output, h5py.h5f.ACC_RDWR, fapl=access_property_list)

    # Reopen the file without a chunk cache
    h5py_file = h5py.File(file_id)

    # Read list of strains
    reader = pd.read_table(kmer_matrix, sep='\t', index_col=0, iterator=True)
    sample_ids = reader.get_chunk(1).columns.values
    logging.debug("The k-mer matrix contains %d samples." % len(sample_ids))

    # Load the strain metadata
    sample_ids, sample_labels = get_mic_dataset(dataset_file=metadata, strains_to_load=sample_ids)
    number_examples = len(sample_ids)

    # Reorder the examples so that the examples of a same class are contiguous
    # This helps count covers faster for the SCM
    reorder = np.argsort(sample_labels)
    sample_ids = sample_ids[reorder]
    sample_labels = sample_labels[reorder]

    # Write strain ids
    logging.debug("Creating the example_identifiers dataset.")
    example_identifiers = h5py_file.create_dataset("example_identifiers",
                                                    data=sample_ids,
                                                    compression=compression,
                                                    compression_opts=compression_opts)

    # Write labels (resistance_labels)
    logging.debug("Creating the labels dataset.")
    labels = h5py_file.create_dataset("labels",
                                    data=sample_labels,
                                    compression=compression,
                                    compression_opts=compression_opts)

    # Initialize kmers (kmer_list) dataset
    logging.debug("Creating the kmers dataset.")
    kmers = h5py_file.create_dataset("kmers",
                                    shape=(kmer_count, ),
                                    dtype=kmer_dtype,
                                    compression=compression,
                                    compression_opts=compression_opts)

    # Initialize attribute_classifications dataset
    logging.debug("Creating the attribute_classifications dataset.")
    attribute_classifications = h5py_file.create_dataset("attribute_classifications",
                                                    shape=(int(ceil(1.0 * number_examples /
                                                                    attribute_classifications_packing_size)),
                                                           kmer_count * 2),
                                                    dtype=attribute_classifications_dtype,
                                                    compression=compression,
                                                    compression_opts=compression_opts,
                                                    chunks=(1, tsv_block_size))  # We will be loading entire lines

    # Initialize type_by_attribute dataset
    logging.debug("Creating the type_by_attribute dataset.")
    type_by_attribute = h5py_file.create_dataset("type_by_attribute",
                                                shape=(kmer_count * 2, ),
                                                dtype=type_by_attribute_dtype,
                                                compression=compression,
                                                compression_opts=compression_opts)

    # Initialize kmer_by_attribute dataset
    logging.debug("Creating the kmer_by_attribute dataset.")
    kmer_by_attribute = h5py_file.create_dataset("kmer_by_attribute",
                                                shape=(kmer_count * 2, ),
                                                dtype=kmer_by_attribute_dtype,
                                                compression=compression,
                                                compression_opts=compression_opts)

    logging.debug("Done creating the datasets.")

    logging.debug("Transferring the data from TSV to the HDF5 datasets.")
    tsv_reader = pd.read_table(kmer_matrix, index_col='kmers', sep='\t', chunksize=tsv_block_size)

    for i, chunk in enumerate(tsv_reader):
        logging.debug("Reading block %d from text file."%(i+1))
        kmers_data = chunk.index.values.astype(kmer_dtype)
        logging.debug("Done reading.")
        read_block_size = kmers_data.shape[0]
        start_block_kmers = i * tsv_block_size                  # i*tsv_block_size
        stop_block_kmers = start_block_kmers + read_block_size   # (i+1)*tsv_block_size or handle the case where the last block is not as long as the others
        start_block_attribute = 2 * start_block_kmers                   # i * 2 * tsv_block_size
        stop_block_attribute = start_block_attribute + 2 * read_block_size   # (i * 2 + 2) * tsv_block_size

        if start_block_kmers > kmer_count:
            break

        # Write kmers
        kmers[start_block_kmers:stop_block_kmers] = kmers_data

        # Write attribute_classifications
        attribute_classification_sorted_by_strains = chunk[sample_ids]
        attribute_classification_data = attribute_classification_sorted_by_strains.T.values.astype(np.uint8)

        attribute_classification_inverse_data = np.zeros(attribute_classification_data.shape, dtype=np.uint8)
        attribute_classification_inverse_data[attribute_classification_data == 0] = 1
        attribute_classification_inverse_data[attribute_classification_data == 1] = 0

        attribute_classifications[:, start_block_attribute:stop_block_attribute] = \
            _pack_binary_bytes_to_ints(np.hstack((attribute_classification_data, attribute_classification_inverse_data)),
                                       pack_size=attribute_classifications_packing_size)
        logging.debug("Done writing block to hdf5.")

        # Write kmer_by_attribute
        kmer_index = np.arange(start_block_kmers, stop_block_kmers, dtype=kmer_by_attribute_dtype)
        kmer_by_attribute[start_block_attribute:stop_block_attribute] = np.hstack((kmer_index, kmer_index))

        # Write type_by_attribute: type=0 for absence and type=1 for presence
        type_by_attribute[start_block_attribute:stop_block_attribute] = \
            np.hstack((np.zeros(read_block_size, dtype=type_by_attribute_dtype),
                       np.ones(read_block_size, dtype=type_by_attribute_dtype)))

        logging.debug("Garbage collection.")
        gc.collect()  # Clear the memory objects created during the iteration, or else the memory will keep growing.

    h5py_file.close()
    logging.debug("Completed.")
