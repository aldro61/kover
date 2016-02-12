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

from functools import partial

from ..utils import _hdf5_open_no_chunk_cache

class KoverDataset(object):
    def __init__(self, file):
        self.dataset_open = partial(_hdf5_open_no_chunk_cache, file)

    @property
    def compression(self):
        dataset = self.dataset_open()
        return dataset.attrs["compression"]

    @property
    def genome_count(self):
        dataset = self.dataset_open()
        return dataset["genome_identifiers"].shape[0]

    @property
    def genome_identifiers(self):
        dataset = self.dataset_open()
        return dataset["genome_identifiers"]

    @property
    def genome_source(self):
        dataset = self.dataset_open()
        return dataset.attrs["genome_source"]

    @property
    def genome_source_type(self):
        dataset = self.dataset_open()
        return dataset.attrs["genome_source_type"]

    @property
    def kmer_by_matrix_column(self):
        dataset = self.dataset_open()
        return dataset["kmer_by_matrix_column"]

    @property
    def kmer_count(self):
        dataset = self.dataset_open()
        return dataset["kmer_sequences"].shape[0]

    @property
    def kmer_length(self):
        dataset = self.dataset_open()
        return len(dataset["kmer_sequences"][0])

    @property
    def kmer_matrix(self):
        dataset = self.dataset_open()
        return dataset["kmer_matrix"]

    @property
    def kmer_sequences(self):
        dataset = self.dataset_open()
        return dataset["kmer_sequences"]

    @property
    def phenotype(self):
        dataset = self.dataset_open()
        return KoverDatasetPhenotype(dataset.attrs["phenotype_name"],
                                     dataset["phenotype"],
                                     dataset.attrs["phenotype_metadata_source"])

    @property
    def splits(self):
        dataset = self.dataset_open()
        return [KoverDatasetSplit(split_name,
                                  split.attrs["train_proportion"],
                                  split["train_genome_idx"],
                                  split["test_genome_idx"],
                                  split["unique_risks"],
                                  split["unique_risk_by_kmer"],
                                  split["unique_risk_by_anti_kmer"],
                                  [KoverDatasetFold(fold_name,
                                                    fold["train_genome_idx"],
                                                    fold["test_genome_idx"],
                                                    fold["unique_risks"],
                                                    fold["unique_risk_by_kmer"],
                                                    fold["unique_risk_by_anti_kmer"]) for fold_name, fold in split["folds"].iteritems()],
                                  split.attrs["random_seed"]) for split_name, split in dataset["splits"].iteritems()]

    @property
    def uuid(self):
        dataset = self.dataset_open()
        return dataset.attrs["uuid"]

    def get_split(self, name):
        dataset = self.dataset_open()
        split = dataset["splits"][name]
        return KoverDatasetSplit(name,
                                 split.attrs["train_proportion"],
                                 split["train_genome_idx"],
                                 split["test_genome_idx"],
                                 split["unique_risks"],
                                 split["unique_risk_by_kmer"],
                                 split["unique_risk_by_anti_kmer"],
                                 [] if not "folds" in split
                                 else [KoverDatasetFold(fold_name,
                                                        fold["train_genome_idx"],
                                                        fold["test_genome_idx"],
                                                        fold["unique_risks"],
                                                        fold["unique_risk_by_kmer"],
                                                        fold["unique_risk_by_anti_kmer"]) for fold_name, fold in split["folds"].iteritems()],
                                 split.attrs["random_seed"])

class KoverDatasetPhenotype(object):
    def __init__(self, name, metadata, metadata_source):
        self.name = name
        self.metadata = metadata
        self.metadata_source = metadata_source

class KoverDatasetSplit(object):
    def __init__(self, name, train_proportion, train_genome_idx, test_genome_idx, unique_risks,
                 unique_risk_by_kmer, unique_risk_by_anti_kmer, folds, random_seed):
        self.name = name
        self.train_proportion = train_proportion
        self.train_genome_idx = train_genome_idx
        self.test_genome_idx = test_genome_idx
        self.unique_risks = unique_risks
        self.unique_risk_by_kmer = unique_risk_by_kmer
        self.unique_risk_by_anti_kmer = unique_risk_by_anti_kmer
        self.folds = folds
        self.random_seed = random_seed

    def __str__(self):
        return "%s   Train genomes: %d (%.3f)   Test genomes: %d (%.3f)   Folds: %d   Random Seed: %d" % \
               (self.name,
                len(self.train_genome_idx),
                self.train_proportion,
                len(self.test_genome_idx),
                1.0 - self.train_proportion,
                len(self.folds),
                self.random_seed)

class KoverDatasetFold(object):
    def __init__(self, name, train_genome_idx, test_genome_idx, unique_risks, unique_risk_by_kmer,
                 unique_risk_by_anti_kmer):
        self.name = name
        self.train_genome_idx = train_genome_idx
        self.test_genome_idx = test_genome_idx
        self.unique_risks = unique_risks
        self.unique_risk_by_kmer = unique_risk_by_kmer
        self.unique_risk_by_anti_kmer = unique_risk_by_anti_kmer
