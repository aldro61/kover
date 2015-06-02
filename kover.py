#!/usr/bin/env python

import argparse
import logging
import sys

class KoverDatasetTool(object):
    def __init__(self):
        self.available_commands = ['convert', 'split']

    def convert(self):
        parser = argparse.ArgumentParser(description='Converts a dataset into the kover dataset format.')
        parser.add_argument('-t', '--type', choices=['tsv'], help='The type of the input file', required=True)
        parser.add_argument('-s', '--samples', help='The k-mer matrix. See the documentation for the appropriate format.', required=True)
        parser.add_argument('-m', '--metadata', help='The metadata file in TSV format. The lines should contain tab separated values and contain the following information: sample_id, class', required=True)
        parser.add_argument('-o', '--output', help='The output file', required=True)
        parser.add_argument('-k', '--kmerlen', type=int, help='The k-mer length. If not specified, will be determined in one pass over the entire dataset (could be long!).')
        parser.add_argument('-c', '--kmercount', type=int, help='The number of k-mers. If not specified, will be determined in one pass over the entire dataset (could be long!).')
        parser.add_argument('-g', '--gzip', type=int, help='The gzip compression level (0 - 9). 0 means no compression. The default value is 4.', default=4)
        parser.add_argument('-v', '--verbose', help='Sets the verbosity level.', default=False, action="store_true")
        args = parser.parse_args(sys.argv[3:])

        if args.verbose:
            logging.basicConfig(level=logging.DEBUG,
                                format="%(asctime)s.%(msecs)d %(levelname)s %(module)s - %(funcName)s: %(message)s")

        if args.type == "tsv":
            from kover.dataset.convert.tsv import convert as convert_tsv
            convert_tsv(kmer_matrix=args.samples,
                        metadata=args.metadata,
                        output=args.output,
                        kmer_len=args.kmerlen,
                        kmer_count=args.kmercount,
                        gzip=args.gzip)

    def split(self):
        parser = argparse.ArgumentParser(description='Splits a kover dataset file into a training set, a testing set and optionally cross-validation folds.')
        parser.add_argument('-d', '--dataset', help='The kover HDF5 dataset to split.', required=True)
        parser.add_argument('-o', '--output', help='The output file.', required=True)
        parser.add_argument('-t', '--trainsize', help='The proportion of the data used for training (default is 0.5).', type=float, default=0.5)
        parser.add_argument('-r', '--randomseed', type=int, help='The random seed used for the split (If not provided a random value between 0 and 4294967295 will be used).')
        parser.add_argument('-f', '--folds', type=int, help='The number of cross-validation folds (default is 0 for none, the minimum value is 2).', default=0)
        parser.add_argument('-g', '--gzip', type=int, help='The gzip compression level (0 - 9). 0 means no compression. The default value is 4.', default=4)
        parser.add_argument('-v', '--verbose', help='Sets the verbosity level.', default=False, action="store_true")
        args = parser.parse_args(sys.argv[3:])

        if args.verbose:
            logging.basicConfig(level=logging.DEBUG,
                                format="%(asctime)s.%(msecs)d %(levelname)s %(module)s - %(funcName)s: %(message)s")

        if args.randomseed is None:
            from random import randint
            args.randomseed = randint(0, 4294967295)

        from kover.dataset.split import split
        split(input=args.dataset,
              output=args.output,
              train_prop=args.trainsize,
              random_seed=args.randomseed,
              n_folds=args.folds,
              gzip=args.gzip)


class Kover(object):

    def __init__(self):
        self.available_commands = ['dataset']

        parser = argparse.ArgumentParser(
            #description='Kover -- Learn from kmers',
            usage='''kover <command> [<args>]

The most commonly used commands are:
    dataset     Dataset manipulation tools
''')
        parser.add_argument('command', help='Command to run')
        # parse_args defaults to [1:] for args, but you need to
        # exclude the rest of the args too, or validation will fail
        args = parser.parse_args(sys.argv[1:2])
        if not hasattr(self, args.command):
            print 'kover: \'%s\' is not a kover command. See \'%s --help\'.' % (args.command, sys.argv[0])

            parser.print_help()
            exit(1)
        # use dispatch pattern to invoke method with same name
        getattr(self, args.command)()

    def dataset(self):
        dataset_tool = KoverDatasetTool()

        parser = argparse.ArgumentParser(description='Dataset manipulation utilities',
                                         usage='''%(prog)s dataset <command> [<args>]
The most commonly used commands are:
    convert     Converts a dataset into the kover dataset format.
    split       Splits a kover dataset file into a training set, a testing set and optionally cross-validation folds.''')

        parser.add_argument('command', help='The dataset manipulation to perform',
                            choices=dataset_tool.available_commands)
        args = parser.parse_args(sys.argv[2:3])
        getattr(dataset_tool, args.command)()


if __name__ == '__main__':
    Kover()