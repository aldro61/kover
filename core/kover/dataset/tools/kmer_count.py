# -*- coding: utf-8 -*-
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

from subprocess import call
from os.path import abspath, dirname, join

def contigs_count_kmers(file_path, out_dir, kmer_size, out_compress, nb_cores, verbose, progress):

	dir_path = dirname(abspath(__file__))

	# Calling MultiDSK tool
	call([str(join(dir_path, "kmer_tools", "multidsk")),
					"-file", str(file_path),
					"-out-dir", str(out_dir),
					"-kmer-size", str(kmer_size),
					"-abundance-min", "1",
					"-out-compress", str(out_compress),
					"-nb-cores", str(nb_cores),
					"-out-tmp", str(out_dir),
					"-verbose", str(verbose),
					"-progress", str(progress)])


def reads_count_kmers(file_path, out_dir, kmer_size, abundance_min, out_compress, nb_cores, verbose, progress):
	dir_path = dirname(abspath(__file__))

	# Calling MultiDSK tool
	call([str(join(dir_path, "kmer_tools", "multidsk")),
					"-file", str(file_path),
					"-out-dir", str(out_dir),
					"-kmer-size", str(kmer_size),
					"-abundance-min", str(abundance_min),
					"-out-compress", str(out_compress),
					"-nb-cores", str(nb_cores),
					"-out-tmp", str(out_dir),
					"-verbose", str(verbose),
					"-progress", str(progress)])
