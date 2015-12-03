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

import numpy as np
cimport numpy as np
cimport cython
from libc.stdint cimport uint32_t, uint64_t


# POPCOUNT 32 bits
cdef extern int __builtin_popcount(unsigned int) nogil

@cython.boundscheck(False)
@cython.wraparound(False)
cdef void _inplace_popcount_32_2d(uint32_t[:, :] arr, uint32_t[:] row_mask) nogil:
    """
    Iterates over the elements of an 2D array and replaces them by their popcount

    Parameters:
    -----------
    arr: numpy_array, dtype=np.uint32, shape=(m, n)
        The array for which the popcounts should be computed.

    row_mask: numpy_array, dtype=np.uint32, shape=(m,)
        A mask used to select the bits to consider in the popcount. The same mask is applied to each row. If a bit is
        set in the mask, the corresponding bit will be considered in the popcount. If a bit is unset in the mask, the
        corresponding bit will not be considered in the popcount.
    """
    cdef int i
    cdef int j
    for i in xrange(arr.shape[0]):
        for j in xrange(arr.shape[1]):
            if arr[i,j] != 0:
                arr[i,j] = __builtin_popcount(arr[i,j] & row_mask[i])

def inplace_popcount_32(arr, row_mask):
    """
    Computes the popcount of each element of a numpy array in-place.
    http://en.wikipedia.org/wiki/Hamming_weight

    Parameters:
    -----------
    arr: numpy_array, dtype=np.uint32, shape=(m, n)
         The array of integers for which the popcounts should be computed.

    row_mask: numpy_array, dtype=np.uint32, shape=(m,)
        A mask used to select the bits to consider in the popcount of each element. The same mask is applied to each
        row. If a bit is set in the mask, the corresponding bit will be considered in the popcount. If a bit is unset
        in the mask, the corresponding bit will not be considered in the popcount.
    """
    _inplace_popcount_32_2d(arr, row_mask)


# POPCOUNT 64 bits

cdef extern int __builtin_popcountl(unsigned long) nogil

@cython.boundscheck(False)
@cython.wraparound(False)
cdef void _inplace_popcount_64_2d(uint64_t[:, :] arr, uint64_t[:] row_mask) nogil:
    """
    Iterates over the elements of an 2D array and replaces them by their popcount

    Parameters:
    -----------
    arr: numpy_array, dtype=np.uint64, shape=(m, n)
        The array for which the popcounts should be computed.

    row_mask: numpy_array, dtype=np.uint64, shape=(m,)
        A mask used to select the bits to consider in the popcount. The same mask is applied to each row. If a bit is
        set in the mask, the corresponding bit will be considered in the popcount. If a bit is unset in the mask, the
        corresponding bit will not be considered in the popcount.
    """
    cdef int i
    cdef int j
    for i in xrange(arr.shape[0]):
        for j in xrange(arr.shape[1]):
            if arr[i,j] != 0:
                arr[i,j] = __builtin_popcountl(arr[i,j] & row_mask[i])

def inplace_popcount_64(arr, row_mask):
    """
    Computes the popcount of each element of a numpy array in-place.
    http://en.wikipedia.org/wiki/Hamming_weight

    Parameters:
    -----------
    arr: numpy_array, dtype=np.uint64, shape=(m, n)
         The array of integers for which the popcounts should be computed.

    row_mask: numpy_array, dtype=np.uint64, shape=(m,)
        A mask used to select the bits to consider in the popcount of each element. The same mask is applied to each
        row. If a bit is set in the mask, the corresponding bit will be considered in the popcount. If a bit is unset
        in the mask, the corresponding bit will not be considered in the popcount.
    """
    _inplace_popcount_64_2d(arr, row_mask)