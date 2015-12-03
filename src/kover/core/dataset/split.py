#!/usr/bin/env python
"""
	Kover: Learn interpretable computational phenotyping models from k-merized genomic data
	Copyright (C) 2015  Alexandre Drouin

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

import h5py as h
import logging
import numpy as np

from h5py.h5f import ACC_RDWR
from math import ceil

from ..utils import _hdf5_open_no_chunk_cache, _minimum_uint_size
from ..learning.set_covering_machine.rules import KmerRuleClassifications

def split(input, identifier, train_prop, random_seed, n_folds, warning_callback=None, error_callback=None,
          progress_callback=None):
    # Execution callback functions
    if warning_callback is None:
        warning_callback = lambda w: logging.warning(w)
    if error_callback is None:
        def normal_raise(exception):
            raise exception
        error_callback = normal_raise
    if progress_callback is None:
        progress_callback = lambda p, m: None

    random_generator = np.random.RandomState(random_seed)

    dataset = _hdf5_open_no_chunk_cache(input, ACC_RDWR)

    if dataset.attrs["phenotype_name"] == "NA":
        raise Exception("A dataset must contain phenotypic metadata to be split.")

    if not "splits" in dataset:
        splits = dataset.create_group("splits")
    else:
        splits = dataset["splits"]

    if identifier in splits:
        raise Exception("A split with the identifier %s already exists in the dataset." % identifier)

    split = splits.create_group(identifier)
    split.attrs["random_seed"] = random_seed
    split.attrs["n_folds"] = n_folds
    split.attrs["train_proportion"] = train_prop

    # for progress
    n_splits_done = 0
    n_splits_to_perform = 1 + n_folds

    # Create the training and testing sets
    logging.debug("Splitting the genomes into a training set and a testing set.")
    n_genomes = len(dataset["genome_identifiers"])
    n_train = int(ceil(train_prop * n_genomes))
    idx = np.arange(n_genomes)
    random_generator.shuffle(idx)
    train_idx = idx[:n_train]
    test_idx = idx[n_train:]
    del idx, n_train
    example_idx_dtype = _minimum_uint_size(n_genomes)
    split.create_dataset("train_genome_idx", data=np.sort(train_idx), dtype=example_idx_dtype)
    split.create_dataset("test_genome_idx", data=np.sort(test_idx), dtype=example_idx_dtype)
    n_splits_done += 0.5
    progress_callback("Split", 1.0 * n_splits_done / n_splits_to_perform)

    # Compute the kmer individual risks (store only a pointer to unique values [rounded at 5 decimals])
    logging.debug("Computing the k-mer individual risks.")
    labels = dataset["phenotype"][...]
    train_pos_idx = train_idx[labels[train_idx] == 1]
    train_neg_idx = train_idx[labels[train_idx] == 0]
    kmer_matrix = KmerRuleClassifications(dataset["kmer_matrix"], len(labels))
    kmer_risks = (len(train_pos_idx) - kmer_matrix.sum_rows(train_pos_idx)[: dataset["kmer_matrix"].shape[1]]).astype(np.float) # n positive errors
    kmer_risks += kmer_matrix.sum_rows(train_neg_idx)[: dataset["kmer_matrix"].shape[1]] # n negative errors
    kmer_risks /= len(train_idx) # n examples
    np.round(kmer_risks, 5, out=kmer_risks)
    anti_kmer_risks = 1.0 - kmer_risks
    np.round(anti_kmer_risks, 5, out=anti_kmer_risks)
    unique_risks, unique_risk_by_kmer_and_antikmer = np.unique(np.hstack((kmer_risks, anti_kmer_risks)), return_inverse=True)
    del kmer_risks, anti_kmer_risks
    split.create_dataset("unique_risks", data=unique_risks)
    split.create_dataset("unique_risk_by_kmer", data=unique_risk_by_kmer_and_antikmer[:dataset["kmer_matrix"].shape[1]], dtype=_minimum_uint_size(len(unique_risks)))
    split.create_dataset("unique_risk_by_anti_kmer", data=unique_risk_by_kmer_and_antikmer[dataset["kmer_matrix"].shape[1]:], dtype=_minimum_uint_size(len(unique_risks)))
    n_splits_done += 0.5
    progress_callback("Split", 1.0 * n_splits_done / n_splits_to_perform)

    if n_folds > 0:
        logging.debug("Splitting the training set into %d cross-validation folds." % n_folds)
        folds = split.create_group("folds")

        # Assign each genome to a fold randomly
        fold_by_training_set_genome = np.arange(len(train_idx)) % n_folds
        np.random.shuffle(fold_by_training_set_genome)

        for fold in xrange(n_folds):
            logging.debug("Fold %d" % (fold + 1))

            fold_group = folds.create_group("fold_%d" % (fold + 1))

            fold_train_idx = train_idx[fold_by_training_set_genome != fold]
            fold_test_idx = train_idx[fold_by_training_set_genome == fold]
            fold_group.create_dataset("train_genome_idx", data=np.sort(fold_train_idx), dtype=example_idx_dtype)
            fold_group.create_dataset("test_genome_idx", data=np.sort(fold_test_idx), dtype=example_idx_dtype)
            n_splits_done += 0.5
            progress_callback("Split", 1.0 * n_splits_done / n_splits_to_perform)

            # Compute the kmer individual risks (store only a pointer to unique values [rounded at 5 decimals])
            logging.debug("Computing the k-mer individual risks.")
            labels = dataset["phenotype"][...]
            fold_train_pos_idx = fold_train_idx[labels[fold_train_idx] == 1]
            fold_train_neg_idx = fold_train_idx[labels[fold_train_idx] == 0]
            kmer_risks = (len(fold_train_pos_idx) - kmer_matrix.sum_rows(fold_train_pos_idx)[: dataset["kmer_matrix"].shape[1]]).astype(np.float) # n positive errors
            kmer_risks += kmer_matrix.sum_rows(fold_train_neg_idx)[: dataset["kmer_matrix"].shape[1]] # n negative errors
            kmer_risks /= len(fold_train_idx) # n examples
            np.round(kmer_risks, 5, out=kmer_risks)
            anti_kmer_risks = 1.0 - kmer_risks
            np.round(anti_kmer_risks, 5, out=anti_kmer_risks)
            unique_risks, unique_risk_by_kmer_and_antikmer = np.unique(np.hstack((kmer_risks, anti_kmer_risks)), return_inverse=True)
            del kmer_risks, anti_kmer_risks
            fold_group.create_dataset("unique_risks", data=unique_risks)
            fold_group.create_dataset("unique_risk_by_kmer", data=unique_risk_by_kmer_and_antikmer[:dataset["kmer_matrix"].shape[1]], dtype=_minimum_uint_size(len(unique_risks)))
            fold_group.create_dataset("unique_risk_by_anti_kmer", data=unique_risk_by_kmer_and_antikmer[dataset["kmer_matrix"].shape[1]:], dtype=_minimum_uint_size(len(unique_risks)))

            n_splits_done += 0.5
            progress_callback("Split", 1.0 * n_splits_done / n_splits_to_perform)
