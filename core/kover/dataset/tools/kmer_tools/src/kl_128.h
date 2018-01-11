/*
*	Kover: Learn interpretable computational phenotyping models from k-merized genomic data
*	Copyright (C) 2015  Alexandre Drouin & Gael Letarte St-Pierre
*
*	This program is free software: you can redistribute it and/or modify
*	it under the terms of the GNU General Public License as published by
*	the Free Software Foundation, either version 3 of the License, or
*	(at your option) any later version.
*
*	This program is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*
*	You should have received a copy of the GNU General Public License
*	along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef KL_128_H
#define KL_128_H

#include <gatb/gatb_core.hpp>
#include <fstream>
#include <hdf5/hdf5.h>
#include <unordered_set>
#include <unordered_map>
#include "callable.h"

const int KMER_MATRIX_PACKING_SIZE_128 = 64;	

class KmerLister128
{
public:
	typedef conditional<KMER_MATRIX_PACKING_SIZE_128==32, unsigned int, unsigned long>::type packing_type;
	
	KmerLister128 (size_t kmerSize);
	void analyse(string input_file, string output_file, string filter, unsigned int compression, unsigned int chunk_size, Callable& progressbar);
	void bit_shift(packing_type* p_buffer, unsigned long* p_nb_kmers);
	string convert(bitset<128> bits);
	
private:
	unordered_map<bitset<128>, unsigned long> k_map;
	unordered_set<bitset<128>> k_set;
	unsigned long nb_kmers = 0;
	unsigned int nb_genomes = 0;
	size_t kmerSize;
	hid_t h5_file;
};
#endif
