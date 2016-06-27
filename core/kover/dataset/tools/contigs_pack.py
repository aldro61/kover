from subprocess import call
from os.path import dirname, abspath

def contigs_pack_kmers(file_path, out_path, filter_singleton, kmer_length, compression, chunk_size, verbose):
	
	dir_path = dirname(abspath(__file__))

	# Calling DSK2Kover tool
	call([dir_path + "/contigs_tools/dsk2kover", 
					"-file", file_path, 
					"-out", out_path, 
					"-filter", filter_singleton, 
					"-kmer-length", kmer_length, 
					"-compression", compression, 
					"-chunk-size", chunk_size,
					"-verbose", verbose])
