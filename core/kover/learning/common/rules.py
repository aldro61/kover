#!/usr/bin/env python
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

import numpy as np

from math import ceil

from .popcount import inplace_popcount_32, inplace_popcount_64
from ...utils import _minimum_uint_size, _unpack_binary_bytes_from_ints

class KmerRule(object):
    def __init__(self, kmer_index, kmer_sequence, type):
        """
        A k-mer rule

        Parameters:
        -----------
        kmer_index: uint
            The index of the k-mer
        kmer_sequence: string
            The nucleotide sequence of the k-mer
        type: string
            The type of rule: presence or absence (use p or a)
        """
        self.kmer_index = kmer_index
        self.kmer_sequence = kmer_sequence
        self.type = type

    def classify(self, X):
        if self.type == "absence":
            return (X[:, self.kmer_index] == 0).astype(np.uint8)
        else:
            return (X[:, self.kmer_index] == 1).astype(np.uint8)

    def inverse(self):
        return KmerRule(kmer_index=self.kmer_index, kmer_sequence=self.kmer_sequence, type="absence" if self.type == "presence" else "presence")

    def __str__(self):
        return ("Absence(" if self.type == "absence" else "Presence(") + self.kmer_sequence + ")"

class LazyKmerRuleList(object):
    """
    By convention, the first half of the list contains presence rules and the second half contains the absence rules in
    the same order.
    """
    def __init__(self, kmer_sequences, kmer_by_rule):
        self.n_rules = kmer_by_rule.shape[0] * 2
        self.kmer_sequences = kmer_sequences
        self.kmer_by_rule = kmer_by_rule

    def __getitem__(self, idx):
        if idx >= self.n_rules:
            raise ValueError("Index %d is out of range for list of size %d" % (idx, self.n_rules))
        if idx >= len(self.kmer_sequences):
            type = "absence"
            kmer_idx = self.kmer_by_rule[idx % len(self.kmer_sequences)]
        else:
            type = "presence"
            kmer_idx = self.kmer_by_rule[idx]
        return KmerRule(idx % len(self.kmer_sequences), self.kmer_sequences[kmer_idx], type)

    def __len__(self):
        return self.n_rules

class BaseRuleClassifications(object):
    def __init__(self):
        pass

    def get_columns(self, columns):
        raise NotImplementedError()

    def remove_rows(self, rows):
        raise NotImplementedError()

    @property
    def shape(self):
        raise NotImplementedError()

    def sum_rows(self, rows):
        raise NotImplementedError()


class KmerRuleClassifications(BaseRuleClassifications):
    """
    Methods involving columns account for presence and absence rules
    """
    # TODO: Clean up. Get rid of the code to handle deleted rows. We don't need this.
    def __init__(self, dataset, n_rows, block_size=None):
        self.dataset = dataset
        self.dataset_initial_n_rows = n_rows
        self.dataset_n_rows = n_rows
        self.dataset_removed_rows = []
        self.dataset_removed_rows_mask = np.zeros(self.dataset_initial_n_rows, dtype=np.bool)
        self.block_size = (None, None)

        if block_size is None:
            if self.dataset.chunks is None:
                self.block_size = (1, self.dataset.shape[1])
            else:
                self.block_size = self.dataset.chunks
        else:
            if len(block_size) != 2 or not isinstance(block_size[0], int) or not isinstance(block_size[1], int):
                raise ValueError("The block size must be a tuple of 2 integers.")
            self.block_size = block_size

        # Get the size of the ints used to store the data
        if self.dataset.dtype == np.uint32:
            self.dataset_pack_size = 32
            self.inplace_popcount = inplace_popcount_32
        elif self.dataset.dtype == np.uint64:
            self.dataset_pack_size = 64
            self.inplace_popcount = inplace_popcount_64
        else:
            raise ValueError("Unsupported data type for packed attribute classifications array. The supported data" +
                             " types are np.uint32 and np.uint64.")

        super(BaseRuleClassifications, self).__init__()

    def get_columns(self, columns):
        """
        Columns can be an integer (or any object that implements __index__) or a sorted list/ndarray.
        """
        #TODO: Support slicing, make this more efficient than getting the columns individually.
        columns_is_int = False
        if hasattr(columns, "__index__"):  # All int types implement the __index__ method (PEP 357)
            columns = [columns.__index__()]
            columns_is_int = True
        elif isinstance(columns, np.ndarray):
            columns = columns.tolist()
        elif isinstance(columns, list):
            pass
        else:
            columns = list(columns)

        # Detect where an inversion is needed (columns corresponding to absence rules)
        columns, invert_result = zip(* (((column if column < self.dataset.shape[1] else column % self.dataset.shape[1]),
                                         (True if column >= self.dataset.shape[1] else False)) for column in columns))
        columns = list(columns)
        invert_result = np.array(invert_result)

        # Don't return rows that have been deleted
        row_mask = np.ones(self.dataset.shape[0] * self.dataset_pack_size, dtype=np.bool)
        row_mask[self.dataset_initial_n_rows:] = False
        row_mask[self.dataset_removed_rows] = False

        # h5py requires that the column indices are sorted
        unique, inverse = np.unique(columns, return_inverse=True)
        result = _unpack_binary_bytes_from_ints(self.dataset[:, unique.tolist()])[row_mask]
        result = result[:, inverse]
        result[:, invert_result] = 1 - result[:, invert_result]

        if columns_is_int:
            return result.reshape(-1)
        else:
            return result

    def remove_rows(self, rows):
        # Find in which dataset the rows must be removed
        dataset_removed_rows = []
        # TODO: This is inefficient!
        for row_idx in rows:
            current_idx = -1
            n_active_elements_seen = 0
            while n_active_elements_seen <= row_idx:
                current_idx += 1
                if not self.dataset_removed_rows_mask[current_idx]:
                    n_active_elements_seen += 1
            dataset_removed_rows.append(current_idx)
        # Update the dataset removed row lists
        # Update the start and stop indexes
        # Adjust the shape
        # Adjust the number of rows in each dataset
        # Store the sorted relative removed row indexes by dataset
        if len(dataset_removed_rows) > 0:
            self.dataset_removed_rows = sorted(set(self.dataset_removed_rows + dataset_removed_rows))
            self.dataset_removed_rows_mask = np.zeros(self.dataset_initial_n_rows, dtype=np.bool)
            self.dataset_removed_rows_mask[self.dataset_removed_rows] = True
            self.dataset_n_rows = self.dataset_initial_n_rows - len(self.dataset_removed_rows)

    @property
    def shape(self):
        return self.dataset_n_rows, self.dataset.shape[1] * 2

    # TODO: allow summing over multiple lists of rows at a time (saves i/o operations)
    def sum_rows(self, rows):
        """
        Note: Assumes that the rows argument does not contain duplicate elements. Rows will not be considered more than once.
        """
        rows = np.asarray(rows)
        result_dtype = _minimum_uint_size(rows.shape[0])
        result = np.zeros(self.dataset.shape[1] * 2, dtype=result_dtype)

        # Builds a mask to turn off the bits of the rows we do not want to count in the sum.
        def build_row_mask(example_idx, n_examples, mask_n_bits):
            if mask_n_bits not in [8, 16, 32, 64, 128]:
                raise ValueError("Unsupported mask format. Use 8, 16, 32, 64 or 128 bits.")

            n_masks = int(ceil(float(n_examples) / mask_n_bits))
            masks = [0] * n_masks

            for idx in example_idx:
                example_mask = idx / mask_n_bits
                example_mask_idx = mask_n_bits - (idx - mask_n_bits * example_mask) - 1
                masks[example_mask] |= 1 << example_mask_idx

            return np.array(masks, dtype="u" + str(mask_n_bits / 8))

        # Find the rows that occur in each dataset and their relative index
        rows = np.sort(rows)
        dataset_relative_rows = []
        for row_idx in rows:
            # Find which row in the dataset corresponds to the requested row
            # TODO: This is inefficient! Could exploit the fact that rows is sorted to reuse previous iterations.
            current_idx = -1
            n_active_elements_seen = 0
            while n_active_elements_seen <= row_idx:
                current_idx += 1
                if not self.dataset_removed_rows_mask[current_idx]:
                    n_active_elements_seen += 1
            dataset_relative_rows.append(current_idx)

        # Create a row mask for each dataset
        row_mask = build_row_mask(dataset_relative_rows, self.dataset_initial_n_rows, self.dataset_pack_size)
        del dataset_relative_rows

        # For each dataset load the rows for which the mask is not 0. Support column slicing aswell
        n_col_blocks = int(ceil(1.0 * self.dataset.shape[1] / self.block_size[1]))
        rows_to_load = np.where(row_mask != 0)[0]
        n_row_blocks = int(ceil(1.0 * len(rows_to_load) / self.block_size[0]))

        for row_block in xrange(n_row_blocks):
            block_row_mask = row_mask[rows_to_load[row_block * self.block_size[0]:(row_block + 1) * self.block_size[0]]]

            for col_block in xrange(n_col_blocks):

                # Load the appropriate rows/columns based on the block sizes
                block = self.dataset[rows_to_load[row_block * self.block_size[0]:(row_block + 1) * self.block_size[0]],
                                     col_block * self.block_size[1]:(col_block + 1) * self.block_size[1]]

                # Popcount
                if len(block.shape) == 1:
                    block = block.reshape(1, -1)
                self.inplace_popcount(block, block_row_mask)

                # Increment the sum
                result[col_block * self.block_size[1]:min((col_block + 1) * self.block_size[1], self.dataset.shape[1])] += np.sum(block, axis=0)

        # Compute the sum for absence rules
        result[self.dataset.shape[1] : ] = len(rows) - result[: self.dataset.shape[1]]

        return result
