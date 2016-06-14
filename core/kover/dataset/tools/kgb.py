from subprocess import call
from os.path import dirname, abspath

def kmers_global_buffering(file_path, out_path, filter_option, kmer_length, compression, chunk_size):
	
	dir_path = dirname(abspath(__file__))

	# Calling KGB tool
	call([dir_path + "/kgb/build/bin/kgb", 
					"-file", file_path, 
					"-out", out_path, 
					"-filter", filter_option, 
					"-kmer-length", kmer_length, 
					"-compression", compression, 
					"-chunk-size", chunk_size])
