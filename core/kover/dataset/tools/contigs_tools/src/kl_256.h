#ifndef KL_256_H
#define KL_256_H

#include <gatb/gatb_core.hpp>
#include <fstream>
#include <hdf5/hdf5.h>
#include <unordered_set>
#include <unordered_map>

const int KMER_MATRIX_PACKING_SIZE_256 = 64;	

class KmerLister256
{
public:
	typedef conditional<KMER_MATRIX_PACKING_SIZE_256==32, unsigned int, unsigned long>::type packing_type;
	
	KmerLister256 (size_t kmerSize);
	void analyse(string input_file, string output_file, string filter, unsigned int compression, unsigned int chunk_size);
	void bit_shift(packing_type* p_buffer, unsigned long* p_nb_kmers);
	string convert(bitset<256> bits);
	
private:
	unordered_map<bitset<256>, unsigned long> k_map;
	unordered_set<bitset<256>> k_set;
	unsigned long nb_kmers = 0;
	unsigned int nb_genomes = 0;
	size_t kmerSize;
	hid_t h5_file;
};
#endif
