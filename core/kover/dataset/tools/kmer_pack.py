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

def contigs_pack_kmers(file_path, out_path, filter_singleton, kmer_length, compression, chunk_size, nb_genomes, progress):

	dir_path = dirname(abspath(__file__))

	# Calling DSK2Kover tool
	call([str(join(dir_path, "kmer_tools", "dsk2kover")),
					"-file", str(file_path),
					"-out", str(out_path),
					"-filter", str(filter_singleton),
					"-kmer-length", str(kmer_length),
					"-compression", str(compression),
					"-chunk-size", str(chunk_size),
					"-nb-genomes", str(nb_genomes),
					"-verbose", str(progress)])


reads_pack_kmers = contigs_pack_kmers
