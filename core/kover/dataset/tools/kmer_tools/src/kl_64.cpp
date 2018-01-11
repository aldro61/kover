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

#include "kl_64.h"

/********************************************************************************/

struct Options_64
{
	Options_64 (Storage* storage, size_t kmerSize,const unordered_map<bitset<64>, unsigned long>& k_map, KmerLister64::packing_type* p_buffer)
	:storage(storage), kmerSize(kmerSize), k_map(k_map), p_buffer(p_buffer) {}
	Storage*   storage;
	size_t     kmerSize;
	const unordered_map<bitset<64>, unsigned long>& k_map;
	KmerLister64::packing_type* p_buffer;
};	

template<size_t span> struct Functor_Analyse_All_64 {void operator () (Options_64 options)
{
	Storage*   storage  = options.storage;
	size_t     kmerSize = options.kmerSize;
	const unordered_map<bitset<64>, unsigned long>& k_map = options.k_map;
	KmerLister64::packing_type* p_buffer = options.p_buffer;
	
	typedef typename Kmer<span>::Count Count;

	// Obtaining solid kmers collection from the 'dsk' group and from the 'solid' collection
	Partition<Count>& solidKmers = storage->getGroup("dsk").getPartition<Count> ("solid");;

	// Iterating over the kmers collection
	Iterator<Count>* itKmers = solidKmers.iterator();
	LOCAL(itKmers);
	for (itKmers->first(); !itKmers->isDone(); itKmers->next())
	{
		const Count& count = itKmers->item();
		
		// Adding presence in the kmers-genome array at the kmer column
		bitset<64> bits (count.value.getVal());
		p_buffer[k_map.at(bits)] += 1;
	}
}};

template<size_t span> struct Functor_Analyse_Filter_64 {void operator () (Options_64 options)
{
	Storage*   storage  = options.storage;
	size_t     kmerSize = options.kmerSize;
	const unordered_map<bitset<64>, unsigned long>& k_map = options.k_map;
	KmerLister64::packing_type* p_buffer = options.p_buffer;
	
	typedef typename Kmer<span>::Count Count;

	// Obtaining solid kmers collection from the 'dsk' group and from the 'solid' collection
	Partition<Count>& solidKmers = storage->getGroup("dsk").getPartition<Count> ("solid");

	// Iterating over the kmers collection
	Iterator<Count>* itKmers = solidKmers.iterator();
	LOCAL(itKmers);
	for (itKmers->first(); !itKmers->isDone(); itKmers->next())
	{
		const Count& count = itKmers->item();
	   
		bitset<64> bits (count.value.getVal());
		// Verifying kmer presence in unordered map to filter singleton
		if (k_map.count(bits))
		{
			// Adding presence in the kmers-genome array at the kmer column
			p_buffer[k_map.at(bits)] += 1;
		}
	}
}};

struct Parameter_64
{
	Parameter_64 (Storage* storage, size_t kmerSize, unordered_map<bitset<64>, unsigned long>& k_map, unsigned long* p_nb_kmers, unordered_set<bitset<64>>& k_set)
	: storage(storage), kmerSize(kmerSize), k_map(k_map), p_nb_kmers(p_nb_kmers), k_set(k_set)  {}
	Storage*   storage;
	size_t     kmerSize;
	unordered_map<bitset<64>, unsigned long>& k_map;
	unsigned long* p_nb_kmers;
	unordered_set<bitset<64>>& k_set;
};

template<size_t span> struct Functor_Read_All_64  {void operator ()  (Parameter_64 parameter)
{
	Storage*   storage  = parameter.storage;
	size_t     kmerSize = parameter.kmerSize;
	unordered_map<bitset<64>, unsigned long>& k_map = parameter.k_map;
	unsigned long* p_nb_kmers = parameter.p_nb_kmers;

	typedef typename Kmer<span>::Count Count;
	
	// Obtaining solid kmers collection from the 'dsk' group and from the 'solid' collection
	Partition<Count>& solidKmers = storage->getGroup("dsk").getPartition<Count> ("solid");
	
	// Iterating over the kmers collection
	Iterator<Count>* itKmers = solidKmers.iterator();
	LOCAL(itKmers);

	for (itKmers->first(); !itKmers->isDone(); itKmers->next())
	{
		const Count& count = itKmers->item();
		// Creating unordered map with every unique kmers
		bitset<64> bits (count.value.getVal());
		if (not(k_map.count(bits)))
		{
			k_map[bits] = *p_nb_kmers;
			(*p_nb_kmers)++;
		}
		
	}
}};

template<size_t span> struct Functor_Read_Filter_64  {void operator ()  (Parameter_64 parameter)
{
	Storage*   storage  = parameter.storage;
	size_t     kmerSize = parameter.kmerSize;
	unordered_map<bitset<64>, unsigned long>& k_map = parameter.k_map;
	unsigned long* p_nb_kmers = parameter.p_nb_kmers;
	unordered_set<bitset<64>>& k_set = parameter.k_set;

	typedef typename Kmer<span>::Count Count;

	// Obtaining solid kmers collection from the 'dsk' group and from the 'solid' collection
	Partition<Count>& solidKmers = storage->getGroup("dsk").getPartition<Count> ("solid");
	
	// Iterating over the kmers collection
	Iterator<Count>* itKmers = solidKmers.iterator();
	LOCAL(itKmers);
	for (itKmers->first(); !itKmers->isDone(); itKmers->next())
	{
		const Count& count = itKmers->item();
		bitset<64> bits (count.value.getVal()); 
		
		/* Creating unordered map while filtering singleton
		 * k_set is used to store kmers that appeared only once
		 * The second time they appear, they are added to the unordered map
		 * */
		if (not(k_map.count(bits)))
		{
			if (k_set.count(bits) == 1)
			{
				k_map[bits] = *p_nb_kmers;
				(*p_nb_kmers)++;
				k_set.erase(bits);
			}
			else
			{ 
				k_set.insert(bits);
			}
		}
	}
}};

/** */
KmerLister64::KmerLister64 (size_t kmerSize):kmerSize(kmerSize)
{}
/** */
void KmerLister64::analyse(string input_file, string output_file, string filter, unsigned int compression, unsigned int chunk_size, Callable& callable)
{
	string line;
	// Opening file containing list of files from dsk to process
	ifstream h5_list (input_file);
	try
	{
		if (h5_list)
		{
			while (getline(h5_list, line))
			{
				nb_genomes ++;
				
				// Loading dsk output file
				// Note : StorageFactory dynamically allocates an instance
				// which explains the use of LOCAL.
				Storage* storage = StorageFactory(STORAGE_HDF5).load (line);
				LOCAL (storage);
				// Launching functor to read kmers with adequate filtering mode
				if (filter == "nothing")
				{
					Integer::apply<Functor_Read_All_64,Parameter_64> (kmerSize, Parameter_64(storage, kmerSize, k_map, &nb_kmers, k_set));
				}
				else if (filter == "singleton")
				{
					Integer::apply<Functor_Read_Filter_64,Parameter_64> (kmerSize, Parameter_64(storage, kmerSize, k_map, &nb_kmers, k_set));	
				}
				callable();
			}
		}
	}
	
	catch (Exception& e)
	{
		h5_list.close();
		throw e;
	}
	// Closing file
	h5_list.close();
	
	// Clearing unordered_set
	k_set.clear();
	
	// Initializing hdf5 interface (dcpl = dataset creation property list)
	hid_t    dataset, datatype, dataspace, memspace, dcpl_matrix;
	
	// Opening output file
	h5_file = H5Fopen((output_file).c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
	hid_t plist;
	plist = H5Fget_access_plist(h5_file);
	H5Pset_cache(plist, 0, 0, 0, 0);
	H5Fclose(h5_file); 
	h5_file = H5Fopen ((output_file).c_str(), H5F_ACC_RDWR, plist);
	
	
	// Initializing compression parameters
	unsigned int chunk_length = (nb_kmers < chunk_size) ? nb_kmers : chunk_size;
	hsize_t chunk[2] = {1, chunk_length};
	
	// Setting dataset creation property list
	dcpl_matrix = H5Pcreate (H5P_DATASET_CREATE);
    H5Pset_deflate (dcpl_matrix, compression);
    H5Pset_chunk (dcpl_matrix, 2, chunk);
	
	// Creating dataspace
	hsize_t dim[2] = {(unsigned long long) (ceil((1.0 * nb_genomes) / KMER_MATRIX_PACKING_SIZE_64)), nb_kmers};
	dataspace = H5Screate_simple(2, dim, NULL); 
	
	// Initializing dataspace subset creation parameters
	hsize_t offset[2] = {0, 0};
	hsize_t stride[2] = {1, 1};
	hsize_t count[2] = {1, nb_kmers};
	hsize_t block[2] = {1, 1};
	
	// Creating memory space for writing
	hsize_t dim_sub[2] = {1, nb_kmers};
	memspace = H5Screate_simple (2, dim_sub, NULL); 

	// Defining datatype
	if (KMER_MATRIX_PACKING_SIZE_64 == 32)
	{
		datatype = H5Tcopy(H5T_NATIVE_UINT);
	}	
	else if (KMER_MATRIX_PACKING_SIZE_64 == 64)
	{
		datatype = H5Tcopy(H5T_NATIVE_ULONG);
	}
	H5Tset_order(datatype, H5T_ORDER_LE);
	
	// Creating dataset in output file
	dataset = H5Dcreate2(h5_file, "kmer_matrix", datatype, dataspace, H5P_DEFAULT, dcpl_matrix, H5P_DEFAULT);
	
	// Creating files name buffer
	string buffer_genomes[KMER_MATRIX_PACKING_SIZE_64];
	
	// Opening file containing list of files from dsk to process
	ifstream dsk_list (input_file);
	try
	{
		if (dsk_list)
		{
			// Creating array data buffer
			packing_type* buffer = new packing_type[nb_kmers];
			for (unsigned int i = 0; i < ceil((1.0 * nb_genomes) / KMER_MATRIX_PACKING_SIZE_64); i++)
			{
				memset(buffer, 0, nb_kmers*sizeof(*buffer));
				unsigned int files_checked = 0;
				bool en_cours = true;
				
				// Filling files name buffer
				while (files_checked < KMER_MATRIX_PACKING_SIZE_64 && en_cours)
				{
					if (getline(dsk_list, line))
					{
						buffer_genomes[files_checked] = line;
						files_checked ++;
					}
					else
					{
						en_cours = false;
					}
				}
				
				// Filling array data buffer
				for (unsigned int genome = 0; genome < files_checked; genome++)
				{
					// Loading dsk output file
					line = buffer_genomes[genome];
					Storage* storage = StorageFactory(STORAGE_HDF5).load (line);
					LOCAL (storage);

					// Launching functor to read kmers and create array with adequate filtering mode
					if (filter == "nothing")
					{
						Integer::apply<Functor_Analyse_All_64,Options_64> (kmerSize, Options_64(storage, kmerSize, k_map, buffer));
					}
					else if (filter == "singleton")
					{
						Integer::apply<Functor_Analyse_Filter_64,Options_64> (kmerSize, Options_64(storage, kmerSize, k_map, buffer));
					}
					callable();
					
					// Shifting array data bits except if last genome
					if (genome != files_checked - 1)
					{
						bit_shift(buffer, &nb_kmers);
					}
					else
					{
						unsigned int empty_rows = KMER_MATRIX_PACKING_SIZE_64 - files_checked;
						if (empty_rows)
						{
							for (unsigned int row = 0; row < empty_rows; row++)
							{
								bit_shift(buffer, &nb_kmers);
							}
						}
					}				
				}
			// Writing array data buffer in output file
			offset[0] = i;
			H5Sselect_hyperslab (dataspace, H5S_SELECT_SET, offset,  stride, count, block);
			H5Dwrite (dataset, datatype, memspace, dataspace, H5P_DEFAULT, buffer);
			}
			delete(buffer);
		}
		// Closing hdf5 interface
		H5Pclose (dcpl_matrix);
		H5Tclose(datatype);
		H5Sclose(memspace);
		H5Sclose(dataspace);
		H5Dclose(dataset);
		
		// Initializing hdf5 interface
		hid_t    dataset_kmers, datatype_kmers, dataspace_kmers, memspace_kmers;
		hid_t 	 dcpl_kmers;
		
		// Defining block size of buffers for output 
		// (number of elements to fill the buffers before a writing operation in the output file)
		unsigned int block_size = 10000;
		
		// Creating buffers for output (k-mer sequences and matrix index)
		char *string_buffer[block_size];
		unsigned long* index = new unsigned long[nb_kmers];
		
		// Initializing dataspace subset creation parameters
		hsize_t offset[] = {0};
		hsize_t stride[] = {1};
		hsize_t count[] = {block_size};
		hsize_t block[] = {1};
		
		// Initializing compression parameters
		hsize_t chunk[] = {block_size};
		
		// Setting dataset creation property list
		dcpl_kmers = H5Pcreate (H5P_DATASET_CREATE);
		H5Pset_deflate (dcpl_kmers, compression);
		H5Pset_chunk (dcpl_kmers, 1, chunk);
		
		// Creating dataspaces
		hsize_t dim[1] = {nb_kmers};
		dataspace_kmers = H5Screate_simple(1, dim, NULL);
		
		// Defining datatypes
		datatype_kmers = H5Tcopy (H5T_C_S1);
		H5Tset_size (datatype_kmers, H5T_VARIABLE);
		
		// Creating memory spaces for writing
		hsize_t dim_sub[1] = {block_size};
		memspace_kmers = H5Screate_simple (1, dim_sub, NULL);  
		
		// Creating datasets in output file
		dataset_kmers = H5Dcreate2(h5_file, "kmer_sequences", datatype_kmers, dataspace_kmers, H5P_DEFAULT, dcpl_kmers, H5P_DEFAULT);
		
		// Iterating over the unordered map
		unordered_map<bitset<64>, unsigned long>::const_iterator iter;
		iter = k_map.begin();
		unsigned int kmer_checked;
		for(kmer_checked = 0; kmer_checked < nb_kmers; kmer_checked++)
		{
			// Position relative to buffer space
			unsigned int element = kmer_checked % block_size;
			
			// Filling buffers
			string_buffer[element] = strdup(convert(iter->first).c_str());
			index[iter->second] = kmer_checked;
			
			// Case buffer is filled
			if ( element == (block_size -1))
			{
				// Writing kmer sequences from buffer in output file
				offset[0] = (kmer_checked / block_size) * block_size;
				H5Sselect_hyperslab (dataspace_kmers, H5S_SELECT_SET, offset,  stride, count, block);
				H5Dwrite (dataset_kmers, datatype_kmers, memspace_kmers, dataspace_kmers, H5P_DEFAULT, string_buffer);
				
				// Releasing string allocated on heap
				for (unsigned int i = 0; i < block_size; i++)
				{
					free(string_buffer[i]);
				}
			}
			iter++;
		}
		
		// Last block of size smaller compare to block_size
		unsigned int last_kmers = nb_kmers % block_size;
		if (not(last_kmers == 0))
		{
			// Adapting subset creation parameters for size
			offset[0] = (kmer_checked / block_size) * block_size;
			count[0] = last_kmers;
			
			// Creating adequate memory spaces
			dim_sub[0] = last_kmers;
			memspace_kmers = H5Screate_simple (1, dim_sub, NULL); 
			
			// Writing kmer sequences from buffer in output file
			H5Sselect_hyperslab (dataspace_kmers, H5S_SELECT_SET, offset,  stride, count, block);
			H5Dwrite (dataset_kmers, datatype_kmers, memspace_kmers, dataspace_kmers, H5P_DEFAULT, string_buffer);
			
			// Releasing string allocated on heap
			for (unsigned int i = 0; i < last_kmers; i++)
			{
				free(string_buffer[i]);
			}
		}
		// Clearing unordered_map
		k_map.clear();
		callable();
		
		// Closing hdf5 interface
		H5Pclose (dcpl_kmers);
		H5Tclose(datatype_kmers);
		H5Sclose(memspace_kmers);
		H5Sclose(dataspace_kmers);
		H5Dclose(dataset_kmers);
		
		hid_t    dataset_index, datatype_index, dataspace_index, memspace_index;
		hid_t 	 dcpl_index;
		
		// Initializing compression parameters
		hsize_t chunk_index[] = {chunk_length};
		
		// Setting dataset creation property list
		dcpl_index = H5Pcreate (H5P_DATASET_CREATE);
		H5Pset_deflate (dcpl_index, compression);
		H5Pset_chunk (dcpl_index, 1, chunk_index);
		
		// Creating dataspaces
		hsize_t dim_index[1] = {nb_kmers};
		dataspace_index = H5Screate_simple(1, dim_index, NULL);
		
		// Defining datatypes
		datatype_index = H5Tcopy(H5T_NATIVE_ULONG);
		H5Tset_order(datatype_index, H5T_ORDER_LE);
		
		// Creating datasets in output file
		dataset_index = H5Dcreate2(h5_file, "kmer_by_matrix_column", datatype_index, dataspace_index, H5P_DEFAULT, dcpl_index, H5P_DEFAULT);
		
		H5Dwrite (dataset_index, datatype_index, H5S_ALL, dataspace_index, H5P_DEFAULT, index);
		
		// Releasing index vector allocated on heap
		delete(index);
		callable();
		
		// Closing hdf5 interface
		H5Pclose (dcpl_index);
		H5Tclose(datatype_index);
		H5Sclose(dataspace_index);
		H5Dclose(dataset_index);
		H5Fclose(h5_file); 
	}
	catch (Exception& e)
	{
		dsk_list.close();
		throw e;
	}
	// Closing file
	dsk_list.close();
}

void KmerLister64::bit_shift(packing_type* p_buffer, unsigned long* p_nb_kmers)
{
	for (unsigned int kmer = 0; kmer < *p_nb_kmers; kmer++)
	{
		p_buffer[kmer] = p_buffer[kmer] << 1;
	}
}

string KmerLister64::convert(bitset<64> bits)
{
	char seq[kmerSize + 1];
	char bin2NT[2][2] = {{'A', 'T'}, {'C', 'G'}};
	
	for (int i = kmerSize - 1; i >= 0; i--)
	{
		seq[i] = bin2NT[bits[0]][bits[1]];
		bits = bits>>2;
	}
	seq[kmerSize] = '\0';
	return seq;
}


