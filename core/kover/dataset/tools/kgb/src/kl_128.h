#ifndef KL_128_H
#define KL_128_H

#include <gatb/gatb_core.hpp>
#include <fstream>
#include <hdf5/hdf5.h>
#include <unordered_set>
#include <unordered_map>

const int KMER_MATRIX_PACKING_SIZE_128 = 32;	

class KmerLister128
{
public:
	typedef conditional<KMER_MATRIX_PACKING_SIZE_128==32, unsigned int, unsigned long>::type packing_type;
	
	KmerLister128 (size_t kmerSize);
	void analyse(string input_file, string output_file, string filter);
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
