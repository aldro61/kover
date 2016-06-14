from subprocess import call
from os.path import dirname, abspath

def kmers_counting_machine(file_path, out_dir, kmer_size, abundance_min, out_compress, nb_cores, verbose):

	dir_path = dirname(abspath(__file__))
	
	# Calling KCM tool
	call([dir_path + "/kgb/build/bin/kcm", 
					"-file", file_path, 
					"-out-dir", out_dir,
					"-kmer-size", kmer_size,
					"-abundance-min", abundance_min,
					"-out-compress", out_compress,
					"-nb-cores", nb_cores,
					"-verbose", verbose])
