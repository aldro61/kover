#!/usr/bin/env python
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

import argparse
import logging
import numpy as np
import warnings; warnings.filterwarnings("ignore")

from collections import defaultdict
from tempfile import gettempdir
from pkg_resources import get_distribution
from sys import argv
from os.path import abspath

KOVER_DESCRIPTION = "Kover: Learn interpretable computational phenotyping models from k-merized genomic data"
VERSION = "2.0.1"

class KoverDatasetCreationTool(object):

    def __init__(self):
        self.available_data_sources = ['from-tsv', 'from-contigs', 'from-reads']

    def from_tsv(self):
        parser = argparse.ArgumentParser(prog="kover dataset create from-tsv",
                                         description='Creates a Kover dataset from genomic data and optionally '
                                                     'phenotypic metadata')
        parser.add_argument('--genomic-data', help='A tab-separated file containing the k-mer matrix.',
                            required=True)
        parser.add_argument('--phenotype-description', help='An informative description that is assigned to the'
                                                            ' phenotypic metadata.')
        parser.add_argument('--phenotype-metadata', help='A file containing the phenotypic metadata.')
        parser.add_argument('--output', help='The Kover dataset to be created.', required=True)
        parser.add_argument('--compression', type=int, help='The gzip compression level (0 - 9). 0 means no compression'
                                                            '. The default value is 4.', default=4)
        parser.add_argument('-x', '--progress', help='Shows a progress bar for the execution.', action='store_true')
        parser.add_argument('-v', '--verbose', help='Sets the verbosity level.', default=False, action='store_true')

        # If no argument has been specified, default to help
        if len(argv) == 4:
            argv.append("--help")

        args = parser.parse_args(argv[4:])

        # Input validation logic
        if (args.phenotype_description is not None and args.phenotype_metadata is None) or (
                        args.phenotype_description is None and args.phenotype_metadata is not None):
            print "Error: The phenotype description and metadata file must be specified."
            exit()

        # Package imports
        from progressbar import Bar, Percentage, ProgressBar, Timer

        if args.verbose:
            logging.basicConfig(level=logging.DEBUG,
                                format="%(asctime)s.%(msecs)d %(levelname)s %(module)s - %(funcName)s: %(message)s")

        if args.progress:
            progress_vars = {"current_task": None, "pbar": None}

            def progress(task_name, p):
                if task_name != progress_vars["current_task"]:
                    if progress_vars["pbar"] is not None:
                        progress_vars["pbar"].finish()
                    progress_vars["current_task"] = task_name
                    progress_vars["pbar"] = ProgressBar(widgets=['%s: ' % task_name, Percentage(), Bar(), Timer()],
                                                        maxval=1.0)
                    progress_vars["pbar"].start()
                else:
                    progress_vars["pbar"].update(p)
        else:
            progress = None

        from kover.dataset.create import from_tsv
        from_tsv(tsv_path=args.genomic_data,
                 output_path=args.output,
                 phenotype_description=args.phenotype_description,
                 phenotype_metadata_path=args.phenotype_metadata,
                 gzip=args.compression,
                 progress_callback=progress)

        if args.progress:
            progress_vars["pbar"].finish()


    def from_contigs(self):
        parser = argparse.ArgumentParser(prog="kover dataset create from-contigs",
                                         description='Creates a Kover dataset from genomic data and optionally '
                                                     'phenotypic metadata')
        parser.add_argument('--genomic-data', help='A tab-separated file with one line per genome in the format '
                                                   'GENOME_ID{tab}PATH, where the path refers to a fasta file '
                                                   'containing the genome\'s contigs.', required=True)
        parser.add_argument('--phenotype-description', help='An informative description that is assigned to the'
                                                            ' phenotypic metadata.')
        parser.add_argument('--phenotype-metadata', help='A file containing the phenotypic metadata.')
        parser.add_argument('--output', help='The Kover dataset to be created.', required=True)
        parser.add_argument('--kmer-size', help='The k-mer size (max is 128). The default is 31.', default=31)
        parser.add_argument('--singleton-kmers', help='Include k-mers that only occur in one genome. Disabled by '
                                                      'default.', default=False, action='store_true')
        parser.add_argument('--n-cpu', '--n-cores', help='The number of cores used by DSK. The default value is 0 (all cores).',
                                                         default=0)
        parser.add_argument('--compression', type=int, help='The gzip compression level (0 - 9). 0 means no compression'
                                                            '. The default value is 4.', default=4)
        parser.add_argument('--temp-dir', help='Output directory for temporary files. The default is the system\'s temp dir.',
                                                default=gettempdir())
        parser.add_argument('-x', '--progress', help='Shows a progress bar for the execution.', action='store_true')
        parser.add_argument('-v', '--verbose', help='Sets the verbosity level.', default=False, action='store_true')

        # If no argument has been specified, default to help
        if len(argv) == 4:
            argv.append("--help")

        args = parser.parse_args(argv[4:])

        # Input validation logic
        if (args.phenotype_description is not None and args.phenotype_metadata is None) or (
                        args.phenotype_description is None and args.phenotype_metadata is not None):
            print "Error: The phenotype description and metadata file must be specified."
            exit()

        if args.verbose:
            logging.basicConfig(level=logging.DEBUG,
                                format="%(asctime)s.%(msecs)d %(levelname)s %(module)s - %(funcName)s: %(message)s")

        from kover.dataset.create import from_contigs

        if not args.singleton_kmers:
            filter_option = "singleton"
        else:
            filter_option = "nothing"

        from_contigs(contig_list_path=args.genomic_data,
                     output_path=args.output,
                     kmer_size=args.kmer_size,
                     filter_singleton=filter_option,
                     phenotype_description=args.phenotype_description,
                     phenotype_metadata_path=args.phenotype_metadata,
                     gzip=args.compression,
                     temp_dir=args.temp_dir,
                     nb_cores=args.n_cpu,
                     verbose=args.verbose,
                     progress=args.progress)

    def from_reads(self):
        parser = argparse.ArgumentParser(prog="kover dataset create from-reads",
                                         description='Creates a Kover dataset from genomic data and optionally '
                                                     'phenotypic metadata')
        parser.add_argument('--genomic-data', help='A tab-separated file with one line per genome in the format '
                                                   'GENOME_ID{tab}PATH, where the path refers to a directory '
                                                   'containing the genome\'s reads in fastq(.gz) files.',
                                                    required=True)
        parser.add_argument('--phenotype-description', help='An informative description that is assigned to the'
                                                            ' phenotypic metadata.')
        parser.add_argument('--phenotype-metadata', help='A file containing the phenotypic metadata.')
        parser.add_argument('--output', help='The Kover dataset to be created.', required=True)
        parser.add_argument('--kmer-size', help='The k-mer size (max is 128). The default is 31.', default=31)
        parser.add_argument('--kmer-min-abundance', help='The minimum number of times a k-mer must be found in a read file '
                                                         'in order to be considered. All k-mers that do not meet this '
                                                         'threshold are discarded. This value should be chosen based on '
                                                         'genome coverage (ex: 100x coverage -> you could use 10). '
                                                         'The default is 1.', default=1)
        parser.add_argument('--singleton-kmers', help='Include k-mers that only occur in one genome. Disabled by '
                                                      'default.', default=False, action='store_true')
        parser.add_argument('--n-cpu', '--n-cores', help='The number of cores used by DSK. The default value is 0 (all cores).',
                                                         default=0)
        parser.add_argument('--compression', type=int, help='The gzip compression level (0 - 9). 0 means no compression'
                                                            '. The default value is 4.', default=4)
        parser.add_argument('--temp-dir', help='Output directory for temporary files. The default is the system\'s temp dir.',
                                                default=gettempdir())
        parser.add_argument('-x', '--progress', help='Shows a progress bar for the execution.', action='store_true')
        parser.add_argument('-v', '--verbose', help='Sets the verbosity level.', default=False, action='store_true')

        # If no argument has been specified, default to help
        if len(argv) == 4:
            argv.append("--help")

        args = parser.parse_args(argv[4:])

        # Input validation logic
        if (args.phenotype_description is not None and args.phenotype_metadata is None) or (
                        args.phenotype_description is None and args.phenotype_metadata is not None):
            print "Error: The phenotype description and metadata file must be specified."
            exit()

        if args.verbose:
            logging.basicConfig(level=logging.DEBUG,
                                format="%(asctime)s.%(msecs)d %(levelname)s %(module)s - %(funcName)s: %(message)s")

        from kover.dataset.create import from_reads

        if not args.singleton_kmers:
            filter_option = "singleton"
        else:
            filter_option = "nothing"

        from_reads(reads_folders_list_path=args.genomic_data,
                   output_path=args.output,
                   kmer_size=args.kmer_size,
                   abundance_min=args.kmer_min_abundance,
                   filter_singleton=filter_option,
                   phenotype_description=args.phenotype_description,
                   phenotype_metadata_path=args.phenotype_metadata,
                   gzip=args.compression,
                   temp_dir=args.temp_dir,
                   nb_cores=args.n_cpu,
                   verbose=args.verbose,
                   progress=args.progress)


class KoverDatasetTool(object):
    def __init__(self):
        self.available_commands = ['create', 'info', 'split']

    def create(self):
        creation_tool = KoverDatasetCreationTool()

        parser = argparse.ArgumentParser(usage= \
'''%(prog)s dataset create <data source> [<args>]
The two available data sources are:
    from-tsv     Create a Kover dataset from genomic data in a tsv format
    from-contigs      Create a Kover dataset from contigs
    from-reads      Create a Kover dataset from reads''')

        parser.add_argument('datasource', help='The type of genomic data to be used.',
                            choices=creation_tool.available_data_sources)

        # If no argument has been specified, default to help
        if len(argv) == 3:
            argv.append("--help")

        args = parser.parse_args(argv[3:4])
        args.datasource = args.datasource.replace("-", "_")
        getattr(creation_tool, args.datasource)()

    def info(self):
        parser = argparse.ArgumentParser(prog="kover dataset info",
                                         description='Prints information about the content of a dataset')
        parser.add_argument('--dataset', help='The Kover dataset for which you require information.', required=True)
        parser.add_argument('--all', help='Prints all the available information.', action='store_true')
        parser.add_argument('--genome-type',
                            help='Prints the type of genomic data that was used to create the dataset.',
                            action='store_true')
        parser.add_argument('--genome-source', help='Prints the source (e.g.: path) from which the genomic data was'
                                                    ' acquired.', action='store_true')
        parser.add_argument('--genome-ids', help='Prints the identifiers of the genomes in the dataset.',
                            action='store_true')
        parser.add_argument('--genome-count', help='Prints the number of genomes in the dataset.', action='store_true')
        parser.add_argument('--kmers', help='Prints the sequence of each k-mer in the dataset (fasta).',
                            action='store_true')
        parser.add_argument('--kmer-len', help='Prints the length of the k-mers in the dataset.', action='store_true')
        parser.add_argument('--kmer-count', help='Prints the number of k-mers in the dataset.', action='store_true')
        parser.add_argument('--phenotype-description', help='Prints the description that was assigned to the phenotype.',
                            action='store_true')
        parser.add_argument('--phenotype-metadata',
                            help='Prints the path of the file from which the phenotypic metadata'
                                 ' was acquired.', action='store_true')
        parser.add_argument('--phenotype-tags', help='Prints the phenotype tags associated to the dataset',
                            action='store_true')
        parser.add_argument('--splits',
                            help='Prints the lists of splits of the dataset that are available for learning.',
                            action='store_true')
        parser.add_argument('--uuid', help='Prints the unique identifier of the dataset.', action='store_true')
        parser.add_argument('--compression', help='Prints the data compression options of the dataset.',
                            action='store_true')
        parser.add_argument('--classification-type', help='Prints the dataset classification type.',
                            action='store_true')

        # If no argument has been specified, default to help
        if len(argv) == 3:
            argv.append("--help")

        args = parser.parse_args(argv[3:])

        # Package imports
        from kover.dataset import KoverDataset

        dataset = args.dataset
        dataset = KoverDataset(dataset)
        #TODO check dataset version for new info types
        if args.genome_type or args.all:
            print "Genome type:", dataset.genome_source_type
            print
        if args.genome_source or args.all:
            print "Genome source:", dataset.genome_source
            print
        if args.genome_ids or args.all:
            print "Genome IDs:"
            for id in dataset.genome_identifiers:
                print id
            print
        if args.genome_count:
            print "Genome count:", dataset.genome_count
            print
        if args.kmers or args.all:
            print "Kmer sequences (fasta):"
            for i, k in enumerate(dataset.kmer_sequences):
                print ">k%d" % (i + 1)
                print k
            print
        if args.kmer_len or args.all:
            print "K-mer length:", dataset.kmer_length
            print
        if args.kmer_count or args.all:
            print "K-mer count:", dataset.kmer_count
            print
        if args.phenotype_description or args.all:
            print "Phenotype description:", dataset.phenotype.description
            print
        if args.phenotype_metadata or args.all:
            if dataset.phenotype.description != "NA":
                print "Phenotype metadata source:", dataset.phenotype.metadata_source
            else:
                print "No phenotype metadata."
            print
        if args.phenotype_tags or args.all:
            print "Phenotype tags: ", ", ".join(dataset.phenotype.tags)
            print
        if args.compression or args.all:
            print "Compression:", dataset.compression
            print
        if args.classification_type or args.all:
            print "Classification type:", dataset.classification_type
            print
        if args.splits or args.all:
            splits = dataset.splits
            if len(splits) > 0:
                print "The following splits are available for learning:"
                for split in splits:
                    print split
            else:
                print "There are no splits available for learning."

    def split(self):
        parser = argparse.ArgumentParser(prog="kover dataset split",
                                         description='Splits a kover dataset file into a training set, a testing set '
                                                     'and optionally cross-validation folds')
        parser.add_argument('--dataset', help='The Kover dataset to be split.', required=True)
        parser.add_argument('--id', help='A unique identifier that will be assigned to the split.', required=True)
        parser.add_argument('--train-size', type=float, help='The proportion of the data that will be reserved for '
                                                             'training the learning algorithm (default is 0.5). '
                                                             ' Alternatively, you can specify which genomes to use '
                                                             'for training and testing by using --train-ids and --test-ids.',
                            default=0.5)
        parser.add_argument('--train-ids', help='File containing the identifiers of the genomes used to train the '
                                                 'learning algorithm. If you provide a value for this argument, you must '
                                                 'also provide a value for --test-ids. File format: one id per line')
        parser.add_argument('--test-ids',  help='File containing the identifiers of the genomes used to evaluate '
                                                'the accuracy of the model generated. If you provide a value for this '
                                                'argument, you must also provide a value for --train-ids. File format: '
                                                'one id per line')
        parser.add_argument('--folds', type=int,
                            help='The number of k-fold cross-validation folds to create (default is 0 for none, '
                                 'the minimum value is 2). Folds are required for using k-fold cross-validation '
                                 'in \'kover learn\'.', default=0)
        parser.add_argument('--random-seed', type=int, help='A random seed used for randomly splitting the data. A '
                                                            'specific seed will always lead to the same split. '
                                                            'If not provided, it is set randomly.')
        parser.add_argument('-v', '--verbose', help='Sets the verbosity level.', default=False, action='store_true')
        parser.add_argument('-x', '--progress', help='Shows a progress bar for the execution.', action='store_true')

        # If no argument has been specified, default to help
        if len(argv) == 3:
            argv.append("--help")

        args = parser.parse_args(argv[3:])

        # Validate the number of cross-validation folds
        if args.folds == 1:
            print "Error: The number of cross-validation folds must be 0 or >= 2."
            exit()

        # Validate that training and testing genome ids are both specified if one of them is specified
        if (args.train_ids is not None and args.test_ids is None) or \
                (args.test_ids is not None and args.train_ids is None):
            print "Error: Training and testing genome identifiers must be specified simultaneously."
            exit()

        # Package imports
        from kover.dataset.split import split_with_ids, split_with_proportion
        from progressbar import Bar, Percentage, ProgressBar, Timer
        from random import randint

        if args.verbose:
            logging.basicConfig(level=logging.DEBUG,
                                format="%(asctime)s.%(msecs)d %(levelname)s %(module)s - %(funcName)s: %(message)s")

        if args.random_seed is None:
            args.random_seed = randint(0, 4294967295)

        if args.progress:
            progress_vars = {"current_task": None, "pbar": None}

            def progress(task_name, p):
                if task_name != progress_vars["current_task"]:
                    if progress_vars["pbar"] is not None:
                        progress_vars["pbar"].finish()
                    progress_vars["current_task"] = task_name
                    progress_vars["pbar"] = ProgressBar(widgets=['%s: ' % task_name, Percentage(), Bar(), Timer()],
                                                        maxval=1.0)
                    progress_vars["pbar"].start()
                else:
                    progress_vars["pbar"].update(p)
        else:
            progress = None

        if args.train_ids is not None and args.test_ids is not None:
            split_with_ids(input=args.dataset,
                           split_name=args.id,
                           train_ids_file=args.train_ids,
                           test_ids_file=args.test_ids,
                           random_seed=args.random_seed,
                           n_folds=args.folds,
                           progress_callback=progress)
        else:
            split_with_proportion(input=args.dataset,
                                  split_name=args.id,
                                  train_prop=args.train_size,
                                  random_seed=args.random_seed,
                                  n_folds=args.folds,
                                  progress_callback=progress)

        if args.progress:
            progress_vars["pbar"].finish()

class KoverLearningTool(object):
    def __init__(self):
            self.available_commands = ['scm', 'tree']

    def scm(self):
        parser = argparse.ArgumentParser(prog='kover learn scm', description='Learn a conjunction/disjunction model using the Set Covering Machine algorithm.')
        parser.add_argument('--dataset', help='The Kover dataset to use for learning.', required=True)
        parser.add_argument('--split', help='The identifier of the split of the dataset to use for learning.',
                            required=True)
        parser.add_argument('--model-type', choices=['conjunction', 'disjunction'], nargs='+',
                            help='Hyperparameter: The type of model to learn, conjunction (logical-AND) or disjunction '
                                 '(logical-OR). You can specify multiple space separated values. Refer to the documentation '
                                 'for more information.',
                            default=['conjunction', 'disjunction'])
        parser.add_argument('--p', type=float, nargs='+',
                            help='Hyperparameter: The value of the trade-off parameter in the rule scoring criterion. '
                                 'You can specify multiple space separated values. Refer to the documentation for '
                                 'more information.',
                            default=[0.1, 0.316, 0.562, 1.0, 1.778, 3.162, 10.0, 999999.0])
        parser.add_argument('--kmer-blacklist', help='A file containing a list of k-mers to remove from the analysis.'
                                                     'These k-mers guaranteed not to be used in the models.'
                                                     'File format: fasta file or text file with one k-mer per line.', required=False)
        parser.add_argument('--max-rules', type=int, help='The maximum number of rules that can be included in the '
                                                          'model.', default=10)
        parser.add_argument('--max-equiv-rules', type=int, help='The maximum number of equivalent rules to report for '
                                                                'each rule in the model. This only affects model interpretation. Use the default unless you'
                                                                ' expect that the rules in the model will be equivalent to more than 10000 other rules.',
                            default=10000)
        parser.add_argument('--hp-choice', choices=['bound', 'cv', 'none'],
                            help='The strategy used to select the best values for hyperparameters. The default is '
                                 'k-fold cross-validation, where k is the number of folds defined in the split. Other '
                                 'strategies, such as bound selection are available (refer to the documentation). Using none '
                                 'selects the first value specified for each hyperparameter.', default='cv')
        parser.add_argument('--bound-max-genome-size', type=int, help='Specify this only if --hp-choice is bound. The '
                                                                      'maximum size, in base pairs, of any genome in the dataset. If you are unsure about this '
                                                                      'value, you should use an overestimation. This will only affect the tightness of the bound '
                                                                      'on the error rate. By default number of k-mers in the dataset is used.')
        parser.add_argument('--random-seed', type=int, help='The random seed used for any random operation. '
                                                            'Set this if only if you require that the same random choices are made between repeats.')
        parser.add_argument('--n-cpu', '--n-cores', type=int, help='The number of CPUs used to select the hyperparameter values. '
                                                      'Make sure your computer has enough RAM to handle multiple simultaneous trainings of the '
                                                      'algorithm and that your storage device will not be a bottleneck (simultaneous reading).',
                            default=1)
        parser.add_argument('--output-dir',
                            help='The directory in which to store Kover\'s output. It will be created if '
                                 'it does not exist.', default='.')
        parser.add_argument('-x', '--progress', help='Shows a progress bar for the execution.', action='store_true')
        parser.add_argument('-v', '--verbose', help='Sets the verbosity level.', default=False, action='store_true')
        parser.add_argument('--authorized-rules', type=str, default="", help=argparse.SUPPRESS)

        # If no argument has been specified, default to help
        if len(argv) == 3:
            argv.append("--help")

        args = parser.parse_args(argv[3:])

        # Package imports
        from datetime import timedelta
        from json import dump as json_dump
        from kover.dataset import KoverDataset
        from kover.learning.experiments.experiment_scm import learn_SCM
        from os import mkdir
        from os.path import abspath, exists, join
        from progressbar import Bar, Percentage, ProgressBar, Timer
        from time import time

        # Input validation
        pre_dataset = KoverDataset(args.dataset)
        dataset_kmer_count = pre_dataset.kmer_count
        classification_type = pre_dataset.classification_type
        phenotype_tags = pre_dataset.phenotype.tags[...]

        # Check if the dataset is compatible
        if classification_type != "binary":
                print "Error: The SCM cannot learn a multi-class classifier"
                exit()

        # - Check that the split exists
        try:
            pre_dataset.get_split(args.split)
        except:
            print "Error: The split (%s) does not exist in the dataset. Use 'kover dataset split' to create it." % args.split
            exit()
        # - Must have at least 2 folds to perform cross-validation
        if args.hp_choice == "cv" and len(pre_dataset.get_split(args.split).folds) < 2:
            print "Error: The split must contain at least 2 folds in order to perform cross-validation. " \
                  "Use 'kover dataset split' to create folds."
            exit()
        del pre_dataset

        if args.verbose:
            logging.basicConfig(level=logging.DEBUG,
                                format="%(asctime)s.%(msecs)d %(levelname)s %(module)s - %(funcName)s [%(process)d]: %(message)s")

        if args.progress:
            progress_vars = {"current_task": None, "pbar": None}

            def progress(task_name, p):
                if task_name != progress_vars["current_task"]:
                    if progress_vars["pbar"] is not None:
                        progress_vars["pbar"].finish()
                    progress_vars["current_task"] = task_name
                    progress_vars["pbar"] = ProgressBar(widgets=['%s: ' % task_name, Percentage(), Bar(), Timer()],
                                                        maxval=1.0)
                    progress_vars["pbar"].start()
                else:
                    progress_vars["pbar"].update(p)
        else:
            progress = None

        args.bound_delta = 0.05  # Fixed value to simplify the user experience

        start_time = time()
        best_hp, best_hp_score, \
        train_metrics, test_metrics, \
        model, rule_importances, \
        equivalent_rules, \
        classifications = learn_SCM(dataset_file=args.dataset,
                                    split_name=args.split,
                                    model_type=args.model_type,
                                    p=args.p,
                                    kmer_blacklist_file=None if args.kmer_blacklist is None else abspath(args.kmer_blacklist),
                                    max_rules=args.max_rules,
                                    max_equiv_rules=args.max_equiv_rules,
                                    bound_delta=args.bound_delta,
                                    bound_max_genome_size=args.bound_max_genome_size
                                    if args.bound_max_genome_size is not None
                                    else dataset_kmer_count,
                                    parameter_selection=args.hp_choice,
                                    n_cpu=args.n_cpu,
                                    random_seed=args.random_seed,
                                    authorized_rules=args.authorized_rules,
                                    progress_callback=progress)
        running_time = timedelta(seconds=time() - start_time)

        if args.progress:
            progress_vars["pbar"].finish()

        # Write a report (user friendly)
        metric_aliases = [("risk", "Error Rate"), ("sensitivity", "Sensitivity"), ("specificity", "Specificity"),
                          ("precision", "Precision"), ("recall", "Recall"), ("f1_score", "F1 Score"),
                          ("tp", "True Positives"), ("tn", "True Negatives"),
                          ("fp", "False Positives"), ("fn", "False Negatives")]
        dataset = KoverDataset(args.dataset)
        split = dataset.get_split(args.split)
        report = ""
        report += "Kover Learning Report\n" + "=" * 21 + "\n"
        report += "\n"
        report += "Running time: %s\n" % running_time
        report += "\n"
        report += "Configuration:\n" + "-" * 14 + "\n"
        for key in sorted(args.__dict__.keys()):
            report += "%s: %s\n" % (str(key), str(args.__dict__[key]))
        report += "\n"
        report += "Data summary:\n" + "-" * 13 + "\n"
        report += "Dataset file: %s\n" % abspath(args.dataset)
        report += "Dataset UUID: %s\n" % dataset.uuid
        report += "Phenotype: %s\n" % dataset.phenotype.description.title()
        report += "Genomic data type: %s\n" % dataset.genome_source_type
        report += "Split: %s\n" % args.split

        report += "Number of genomes used for training: %d " % (len(split.train_genome_idx))
        nb_genome_training = {c:(dataset.phenotype.metadata[split.train_genome_idx] == c).sum()\
                                for c in range(len(phenotype_tags))}
        training_groups = ["Group %s: %d" % (phenotype_tags[c], nb_genome_training[c]) for c in range(len(phenotype_tags))]
        report += "(%s)\n" % ", ".join(training_groups)

        report += "Number of genomes used for testing: %d " % (len(split.test_genome_idx))
        nb_genome_testing = {c:(dataset.phenotype.metadata[split.test_genome_idx] == c).sum() if \
                            len(split.test_genome_idx) > 0 else 0 for c in range(len(phenotype_tags))}
        testing_groups = ["Group %s: %d" % (phenotype_tags[c], nb_genome_testing[c]) for c in range(len(phenotype_tags))]
        report += "(%s)\n" % ", ".join(testing_groups)

        report += "Number of k-mers: %d\n" % dataset.kmer_count
        if dataset.genome_source_type == "contigs":
            report += "K-mer size : %s\n" % dataset.kmer_length
            report += "K-mer filtering : %s\n" % dataset.kmer_filter
        report += "\n"
        report += "Hyperparameter Values:\n" + "-" * 22 + "\n"
        if args.hp_choice == "cv":
            report += "Selection strategy: %d-fold cross-validation (score = %.5f)\n" % (
                len(split.folds), best_hp_score)
        elif args.hp_choice == "bound":
            report += "Selection strategy: bound selection (score = %.5f)\n" % best_hp_score
        else:
            report += "Selection strategy: No selection\n"
        report += "Model type: %s\n" % best_hp["model_type"]
        report += "p: %f\n" % best_hp["p"]
        report += "Maximum number of rules: %d\n" % best_hp["max_rules"]
        report += "\n"
        # Write the training set metrics
        report += "Metrics (training data)\n" + "-" * 23 + "\n"
        for key, alias in metric_aliases:
            report += "%s: %s\n" % (str(alias), str(round(train_metrics[key][0], 5)))
        report += "\n"
        # Write the testing set metrics
        if test_metrics is not None:
            report += "Metrics (testing data)\n" + "-" * 22 + "\n"
            for key, alias in metric_aliases:
                report += "%s: %s\n" % (str(alias), str(round(test_metrics[key][0], 5)))
            report += "\n"
        # Write the model
        report += "Model (%s - %d rules):\n" % (model.type.title(), len(model)) + "-" * (
            18 + len(model.type) + len(str(len(model)))) + "\n"
        report += ("\n%s\n" % ("AND" if model.type == "conjunction" else "OR")).join(
            ["%s [Importance: %.2f, %d equivalent rules]" % (str(rule), importance, len(equivalent_rules[i])) for
             i, (rule, importance) in enumerate(zip(model, rule_importances))])
        report += "\n"

        # Show the report
        print report

        # Create the output directory
        if not exists(args.output_dir):
            mkdir(args.output_dir)

        # Write report
        with open(join(args.output_dir, "report.txt"), "w") as f:
            f.write(report)

        # Save detailed results to json
        results = {"data": {"uuid": dataset.uuid,
                            "path": dataset.path,
                            "split": args.split},
                   "cv": {"best_hp": {"values": best_hp,
                                      "score": best_hp_score},
                          "candidate_hp": {"model_type": args.model_type,
                                           "p": args.p,
                                           "max_rules": args.max_rules},
                          "strategy": args.hp_choice},
                   "metrics": {"train": train_metrics,
                               "test": test_metrics},
                   "model": {"n_rules": len(model),
                             "rules": [str(r) for r in model],
                             "rule_importances": rule_importances.tolist(),
                             "type": best_hp["model_type"]},
                   "classifications": classifications,
                   "running_time": running_time.seconds}
        with open(join(args.output_dir, 'results.json'), 'w') as f:
            json_dump(results, f)

        # Save command line arguments to json
        config = dict((k, v) for k, v in args.__dict__.iteritems())
        config["dataset"] = abspath(config['dataset'])
        with open(join(args.output_dir, 'config.json'), 'w') as f:
            json_dump(config, f)

        # Save model (also equivalent rules) [Fasta]
        with open(join(args.output_dir, 'model.fasta'), "w") as f:
            for i, (rule, importance) in enumerate(zip(model, rule_importances)):
                f.write(">rule-%d %s, importance: %.2f\n%s\n\n" % (i + 1, rule.type, importance, rule.kmer_sequence))
                with open(join(args.output_dir, "model_rule_%i_equiv.fasta" % (i + 1)), "w") as f_equiv:
                    f_equiv.write("\n\n".join(
                        [">rule-%d-equiv-%d,%s\n%s" % (i + 1, j + 1, rule.type, rule.kmer_sequence)
                         for j, rule in enumerate(equivalent_rules[i])]))

    def tree(self):
        parser = argparse.ArgumentParser(prog='kover learn tree', description='Learn a decision tree model using the Classification And Regression Trees algorithm.')
        parser.add_argument('--dataset', help='The Kover dataset to use for learning.', required=True)
        parser.add_argument('--split', help='The identifier of the split of the dataset to use for learning.',
                            required=True)
        parser.add_argument('--criterion', type=str, nargs='+', help='Hyperparameter: The criterion used to partition the leaves of the decision tree. '
                            'You can specify multiple space separated values. Refer to the documentation '
                            'for more information. (default: gini)',
                            choices=['gini', 'crossentropy'], default='gini', required=False)
        parser.add_argument('--max-depth', type=int, nargs='+', help='Hyperparameter: The maximum depth of the decision tree. You can specify multiple space separated '
                            'values. Refer to the documentation for more information. (default: 10)',
                            default=10, required=False)
        parser.add_argument('--min-samples-split', type=int, nargs='+',
                            help='Hyperparameter: The minimum number of genomes that a leaf must contain to be partionned into two new leaves. You can specify multiple '
                            'space separated values. Refer to the documentation for more information. (default: 2)', default=2, required=False)
        parser.add_argument('--class-importance', type=str, nargs='+', help='Hyperparameter (this one is tricky so read carefully): This controls the cost of making prediction '
                            'errors on each class. You can either specify a set of space-separated values, which will result in trying each value for each class, or provide a '
                            'set of class-specific values to try using the following syntax: "class1: v1 v2 class2: v1 v2 v3 ...". Refer to the documentation for more information. '
                            '(default: 1.0 for each class)',
                            default=None, required=False)
        parser.add_argument('--kmer-blacklist', help='A file containing a list of k-mers to remove from the analysis.'
                                                     'These k-mers guaranteed not to be used in the models.'
                                                     'File format: fasta file or text file with one k-mer per line.', required=False)
        parser.add_argument('--hp-choice', choices=['bound', 'cv'],
                            help='The strategy used to select the best values for the hyperparameters. The default is '
                                 'k-fold cross-validation, where k is the number of folds defined in the split. '
                                 'Alternatively, you can use bound selection, which is much faster, but may not work '
                                 'as well for smaller datasets. Also note that, even if a single value is provided '
                                 'for each hyperparameter, a search will be performed as part of the tree '
                                 'pruning process. Refer to the documentation for more information. (default: cv)', default='cv')
        parser.add_argument('--bound-max-genome-size', type=int, help='Specify this only if --hp-choice is bound. The '
                            'maximum size, in base pairs, of any genome in the dataset. If you are unsure about this '
                            'value, you should use an overestimation. This will only affect the tightness of the bound '
                            'on the error rate. By default number of k-mers in the dataset is used.')
        parser.add_argument('--n-cpu', '--n-cores', type=int, help='The number of CPUs used to select the hyperparameter values. '
                                                      'Make sure your computer has enough RAM to handle multiple simultaneous trainings of the '
                                                      'algorithm and that your storage device will not be a bottleneck (simultaneous reading).',
                            default=1)
        parser.add_argument('--output-dir',
                            help='The directory in which to store Kover\'s output. It will be created if '
                                 'it does not exist.', default='.')
        parser.add_argument('-x', '--progress', help='Shows a progress bar for the execution.', action='store_true')
        parser.add_argument('-v', '--verbose', help='Sets the verbosity level.', default=False, action='store_true'),
        parser.add_argument('--authorized-rules', type=str, default="", help=argparse.SUPPRESS)  # Hidden argument

        # If no argument has been specified, default to help
        if len(argv) == 3:
            argv.append("--help")

        args = parser.parse_args(argv[3:])

        # Package imports
        from datetime import timedelta
        from json import dump as json_dump
        from kover.dataset import KoverDataset
        from kover.learning.experiments.experiment_cart import learn_CART
        from os import mkdir
        from os.path import abspath, exists, join
        from progressbar import Bar, Percentage, ProgressBar, Timer
        from time import time
        from itertools import permutations, product

        # Input validation
        pre_dataset = KoverDataset(args.dataset)

        # - Check that the split exists
        try:
            pre_dataset.get_split(args.split)
        except:
            print "Error: The split (%s) does not exist in the dataset. Use 'kover dataset split' to create it." % args.split
            exit()


        # - Must have at least 2 folds to perform cross-validation
        if args.hp_choice == "cv" and len(pre_dataset.get_split(args.split).folds) < 2:
            print "Error: The split must contain at least 2 folds in order to perform cross-validation. " \
                  "Use 'kover dataset split' to create folds."
            exit()

        # Load classification task specifications
        dataset_kmer_count = pre_dataset.kmer_matrix.shape[1]
        phenotype_tags = pre_dataset.phenotype.tags[...].tolist()
        classification_type = pre_dataset.classification_type
        del pre_dataset

        def parse_class_importances(class_importance_input):
            """
            Parses a class importance specification with the following format:

            class_1: v11 v12 class_2: v21 v22 ...

            """
            # ---- Input validation follows ----

            # Check that importance values are specified for each class
            for class_name in phenotype_tags:
                try:
                    class_importance_input.index(class_name + ":")
                except ValueError:
                    print "Error: no class importances defined for class \"{}\" which is in the dataset.".format(class_name)
                    exit()

            # Check that no unknown classes were used
            for class_name in [x.replace(":", "") for x in class_importance_input if x.endswith(":")]:
                try:
                    phenotype_tags.index(class_name)
                except ValueError:
                    print "Error: unknown class \"{}\" in class importances.".format(class_name)
                    exit()

            # Check that there is at least one importance value for each class
            for i in xrange(len(class_importance_input)):
                if class_importance_input[i].endswith(":"):
                    try:
                        assert not class_importance_input[i + 1].endswith(":")
                    except:
                        print "Error: no class importances defined for class {} which is in the dataset.".format(class_importance_input[i].replace(":", ""))
                        exit()

            # Check that all importances specified are valid floats
            for v in class_importance_input:
                if not v.endswith(":"):
                    try:
                        float(v)
                    except ValueError:
                        print "Error: invalid value \"{}\" encountered in class importances.".format(v)
                        exit()

            # ---- Input is valid ----

            # Get the list of importances to try for each class
            class_importances = defaultdict(list)
            for v in class_importance_input:
                if v.endswith(":"):
                    current_class = v.replace(":", "")
                else:
                    class_importances[phenotype_tags.index(current_class)].append(float(v))

            # Make a list of dictionnaries giving every possible combination of importances
            grid_classes = class_importances.keys()
            grid = product(*class_importances.values())
            class_importances = [{c: importance for c, importance in zip(grid_classes, grid_row)} for grid_row in grid]

            return np.unique(class_importances)

        # Determine which type of class importances was received as input
        if args.class_importance:
            if args.class_importance[0].endswith(":"):
                # Separate set of class importances for each class
                # Format: c1: v11 v2 c12: v21 v22 ...
                class_importances = parse_class_importances(args.class_importance)
            else:
                # Identical set of class importances for each class
                # Format: v1 v2 v3 ...
                tmp = []
                for c in phenotype_tags:
                    tmp.append(c + ":")
                    tmp += args.class_importance
                class_importances = parse_class_importances(tmp)
        else:
            # No class importances specified, so each has an importance of 1.0
            class_importances = [{c: 1.0 for c in range(len(phenotype_tags))}]

        if args.verbose:
            logging.basicConfig(level=logging.DEBUG,
                                format="%(asctime)s.%(msecs)d %(levelname)s %(module)s - %(funcName)s [%(process)d]: %(message)s")

        if args.progress:
            progress_vars = {"current_task": None, "pbar": None}
            def progress(task_name, p):
                if task_name != progress_vars["current_task"]:
                    if progress_vars["pbar"] is not None:
                        progress_vars["pbar"].finish()
                    progress_vars["current_task"] = task_name
                    progress_vars["pbar"] = ProgressBar(widgets=['%s: ' % task_name, Percentage(), Bar(), Timer()], maxval=1.0)
                    progress_vars["pbar"].start()
                else:
                    progress_vars["pbar"].update(p)
        else:
            progress = None

        args.bound_delta = 0.05  # Fixed value to simplify the user experience
        if args.bound_max_genome_size is None:
            args.bound_max_genome_size = dataset_kmer_count

        start_time = time()
        best_hp, best_hp_score, \
        train_metrics, test_metrics, \
        model, rule_importances, \
        equivalent_rules, \
        classifications = learn_CART(dataset_file=args.dataset,
                                split_name=args.split,
                                criterion=args.criterion,
                                max_depth=args.max_depth,
                                min_samples_split=args.min_samples_split,
                                class_importance=class_importances,
                                bound_delta=args.bound_delta,
                                bound_max_genome_size=args.bound_max_genome_size,
                                kmer_blacklist_file=None if args.kmer_blacklist is None else abspath(args.kmer_blacklist),
                                parameter_selection=args.hp_choice,
                                authorized_rules=args.authorized_rules,
                                n_cpu=args.n_cpu,
                                progress_callback=progress)
        running_time = timedelta(seconds=time() - start_time)

        if args.progress:
            progress_vars["pbar"].finish()

        # Write a report (user friendly)
        if classification_type == "binary":
            metric_aliases = [("risk", "Error Rate"), ("sensitivity", "Sensitivity"), ("specificity", "Specificity"),
                              ("precision", "Precision"), ("recall", "Recall"), ("f1_score", "F1 Score"),
                              ("tp", "True Positives"), ("tn", "True Negatives"),
                              ("fp", "False Positives"), ("fn", "False Negatives")]
        elif classification_type == "multiclass":
            metric_aliases = [("risk", "Error rate"), ("confusion_matrix", "Confusion Matrix")]

        # Convert confusion matrix to a nice text output
        def confusion_matrix_to_str(confusion_matrix):
            matrix_str = ""
            size_header = len(max(phenotype_tags, key=len)) + 5
            col_width = 5
            horizontal_bar = "+-" + "-" * size_header + "+" + "+".join(["-" * col_width for c in range(len(phenotype_tags))]) + "+\n"

            matrix_str += horizontal_bar
            matrix_str += "| " + " " * size_header + "|"
            matrix_str += "|".join([str(c).center(col_width) for c in range(len(phenotype_tags))])
            matrix_str += "|\n"
            matrix_str += horizontal_bar.replace('-', '=')
            for c in range(len(phenotype_tags)):
                matrix_str += "| " + phenotype_tags[c].ljust(size_header - 5) + ("(" + str(c) + ")").center(5) + "|"
                matrix_str += "|".join([str(score).center(col_width) for score in confusion_matrix[c]])  + "|\n"
                matrix_str += horizontal_bar

            return matrix_str

        # Create identifiers for each rules
        rule_ids = {}
        id_by_node = {n: "{}___ex_{}___eq_{}".format(i, n.n_examples, len(equivalent_rules[n.rule])) if not n.is_leaf else "leaf___ex_{}___{}".format(n.n_examples, "__".join(["{0!s}_{1:d}_{2:.8f}".format(model.class_tags[c], len(n.class_examples_idx[c]), n.breiman_info.p_j_given_t[c]) for c in n.class_proportions])) for i, n in model.decision_tree}
        for node_id, node in model.decision_tree:
            if not node.is_leaf:
                simple_rule_id = str(node_id)
                fasta_header = "rule_id: {}, left_child: {}, right_child: {}".format(id_by_node[node], id_by_node[node.left_child], id_by_node[node.right_child])
                rule_ids[node.rule] = {"simple": simple_rule_id, "fasta": fasta_header}
        del id_by_node

        dataset = KoverDataset(args.dataset)
        split = dataset.get_split(args.split)

        report = ""
        report += "Kover Learning Report\n" + "="*21 + "\n"
        report += "\n"
        report += "Running time: %s\n" % running_time
        report += "\n"
        report += "Configuration:\n" + "-"*14 + "\n"
        for key in sorted(args.__dict__.keys()):
            report +=  "%s: %s\n" % (str(key), str(args.__dict__[key]))
        report += "\n"
        report += "Data summary:\n" + "-"*13 + "\n"
        report += "Dataset file: %s\n" % abspath(args.dataset)
        report += "Dataset UUID: %s\n" % dataset.uuid
        report += "Phenotype: %s\n" % dataset.phenotype.description.title()
        report += "Split: %s\n" % args.split

        report += "Number of genomes used for training: %d " % (len(split.train_genome_idx))
        nb_genome_training = {c:(dataset.phenotype.metadata[split.train_genome_idx] == c).sum()\
                                for c in range(len(phenotype_tags))}
        training_groups = ["Group %s: %d" % (phenotype_tags[c], nb_genome_training[c]) for c in range(len(phenotype_tags))]
        report += "(%s)\n" % ", ".join(training_groups)

        report += "Number of genomes used for testing: %d " % (len(split.test_genome_idx))
        nb_genome_testing = {c:(dataset.phenotype.metadata[split.test_genome_idx] == c).sum() if \
                            len(split.test_genome_idx) > 0 else 0 for c in range(len(phenotype_tags))}
        testing_groups = ["Group %s: %d" % (phenotype_tags[c], nb_genome_testing[c]) for c in range(len(phenotype_tags))]
        report += "(%s)\n" % ", ".join(testing_groups)

        report += "\n"
        report += "Hyperparameter Values:\n" + "-" * 22 + "\n"
        if args.hp_choice == "cv":
            report += "Selection strategy: {0:d}-fold cross-validation (score = {1:.5f})\n".format(len(split.folds), best_hp_score)
        elif args.hp_choice == "bound":
            report += "Selection strategy: sample-compression bound (delta = {0:.3f}, max-genome-size = {1:d}, value = {2:.5f})\n".format(args.bound_delta, args.bound_max_genome_size, best_hp_score)
        else:
            report += "Selection strategy: No selection\n"
        report += "Criterion: {0!s}\n".format(best_hp["criterion"])
        report += "Class importance: {0!s}\n".format(", ".join(["class {0!s}: {1:.3f}".format(c, best_hp["class_importance"][i]) for i, c in enumerate(phenotype_tags)]))
        report += "Maximum tree depth: {0:d}\n".format(best_hp["max_depth"])
        report += "Minimum samples to split a node (examples): {0:.3f}\n".format(best_hp["min_samples_split"])
        report += "Pruning alpha: {0:.8f}\n".format(best_hp["pruning_alpha"])
        report += "\n"
        # Print the training set metrics
        report += "Metrics (training data)\n" + "-"*23 + "\n"
        for key, alias in metric_aliases:
            if key == "confusion_matrix":
                report += "%s :\n%s\n" % (str(alias), confusion_matrix_to_str(train_metrics[key][0]))
            else:
                report += "%s: %s\n" % (str(alias), str(round(train_metrics[key][0], 5)))
        report += "\n"
        # Print the testing set metrics
        if test_metrics is not None:
            report += "Metrics (testing data)\n" + "-"*22 + "\n"
            for key, alias in metric_aliases:
                if key == "confusion_matrix":
                    report += "%s :\n%s\n" % (str(alias), confusion_matrix_to_str(test_metrics[key][0]))
                else:
                    report += "%s: %s\n" % (str(alias), str(round(test_metrics[key][0], 5)))
            report += "\n"
        report += "Model (%d rules, depth = %d):\n" % (len(model.decision_tree.rules), model.depth)
        report += str(model) + "\n"
        report += "\n"

        print report

        if not exists(args.output_dir):
            mkdir(args.output_dir)

        # Write report
        with open(join(args.output_dir, "report.txt"), "w") as f:
            f.write(report)

        # Write metrics (not user friendly) [json]
        results = {"data": {"uuid": dataset.uuid,
                            "path": dataset.path,
                            "split": args.split},
                   "cv": {"best_hp": {"values": best_hp,
                                      "score": best_hp_score},
                          "candidate_hp": {"criterion": args.criterion,
                                           "max_depth": args.max_depth},
                          "strategy": args.hp_choice},
                   "metrics": {"train": train_metrics,
                               "test": test_metrics},
                   "model": {"n_rules": len(model.decision_tree.rules),
                             "depth": model.depth,
                             "rules": [str(r) for r in model.decision_tree.rules],
                             "rule_importances": [rule_importances[r] for r in model.decision_tree.rules],
                             "rule_identifiers": [rule_ids[r]["simple"] for r in model.decision_tree.rules]},
                   "classifications": classifications,
                   "running_time": running_time.seconds}
        with open(join(args.output_dir, 'results.json'), 'w') as f:
            json_dump(results, f)

        # Save command line arguments to json
        config = dict(args.__dict__)
        config["dataset"] = abspath(config['dataset'])
        with open(join(args.output_dir, 'config.json'), 'w') as f:
            json_dump(config, f)

        # Save model (also equivalent rules) [json]
        with open(join(args.output_dir, 'model.fasta'), "w") as f:
            for i, rule in enumerate(model.decision_tree.rules):
                f.write(">{0!s}, importance: {1:.2f}\n{2!s}\n\n".format(rule_ids[rule]["fasta"], rule_importances[rule], rule.kmer_sequence))

                with open(join(args.output_dir, "model_rule_{0!s}_equiv.fasta".format(rule_ids[rule]["simple"])), "w") as f_equiv:
                    f_equiv.write("\n\n".join(
                        [">rule-{0!s}-equiv-{1:d}\n{2!s}".format(rule_ids[rule]["simple"], j + 1, equiv_rule.kmer_sequence)
                         for j, equiv_rule in enumerate(equivalent_rules[rule])]))


class CommandLineInterface(object):
    def __init__(self):
        self.available_commands = ['dataset', 'learn']

        parser = argparse.ArgumentParser(description=KOVER_DESCRIPTION)
        parser.add_argument('--cite', action='store_true',
                            help='Show the references to cite if you use Kover in your work')
        parser.add_argument('--license', action='store_true',
                            help='Show the license for your release of Kover')
        parser.add_argument('--version', action='store_true',
                            help='Show the version of your release of Kover')
        parser.add_argument('command', help='Available commands', choices=self.available_commands)

        # If no argument has been specified, default to help
        if len(argv) == 1:
            argv.append("--help")

        # Detect if an option or a command was specified
        if "--" == argv[1][:2]:
            args = parser.parse_args([argv[1], "learn"])

            if args.license:
                print \
"""
Kover: Learn interpretable computational phenotyping models from k-merized genomic data
Copyright (C) 2019  Alexandre Drouin

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
            elif args.version:
                print "cli-%s" % VERSION
                print "core-%s" % get_distribution('kover').version

            elif args.cite:
                print \
"""If you used Kover in your work, please cite our papers:

Drouin, A., Letarte, G., Raymond, F., Marchand, M., Corbeil, J., & Laviolette, F. (2019). Interpretable genotype-to-phenotype classifiers with performance guarantees. Scientific reports, 9(1), 4071.

Drouin, A., Giguère, S., Déraspe, M., Marchand, M., Tyers, M., Loo, V. G., Bourgault, A. M., Laviolette, F. & Corbeil, J. (2016). Predictive computational phenotyping and biomarker discovery using reference-free genome comparisons. BMC Genomics, 17(1), 754.


Bibtex:

@article{Drouin2019,
  title={Interpretable genotype-to-phenotype classifiers with performance guarantees},
  author={Drouin, Alexandre and Letarte, Ga{\"e}l and Raymond, Fr{\'e}d{\'e}ric and Marchand, Mario and Corbeil, Jacques and Laviolette, Fran{\c{c}}ois},
  journal={Scientific reports},
  volume={9},
  number={1},
  pages={4071},
  year={2019},
  publisher={Nature Publishing Group}
}

@article{Drouin2016,
author="Drouin, Alexandre
and Gigu{\\`e}re, S{\\'e}bastien
and D{\\'e}raspe, Maxime
and Marchand, Mario
and Tyers, Michael
and Loo, Vivian G.
and Bourgault, Anne-Marie
and Laviolette, Fran{\\c{c}}ois
and Corbeil, Jacques",
title="Predictive computational phenotyping and biomarker discovery using reference-free genome comparisons",
journal="BMC Genomics",
year="2016",
volume="17",
number="1",
pages="754",
issn="1471-2164",
doi="10.1186/s12864-016-2889-6",
url="http://dx.doi.org/10.1186/s12864-016-2889-6"
}
"""
        else:
            args = parser.parse_args(argv[1:2])
            if not hasattr(self, args.command):
                print 'kover: \'%s\' is not a kover command. See \'%s --help\'.' % (args.command, argv[0])

                parser.print_help()
                exit(1)
            getattr(self, args.command)()

    def dataset(self):
        dataset_tool = KoverDatasetTool()

        parser = argparse.ArgumentParser(usage= \
'''%(prog)s dataset <command> [<args>]
The most commonly used commands are:
    create     Create Kover datasets from genomic data
    split      Split a Kover dataset file into a training set, a testing set and optionally cross-validation folds
    info       Get information about the content of a Kover dataset''')

        parser.add_argument('command', help='The dataset manipulation to perform',
                            choices=dataset_tool.available_commands)

        # If no argument has been specified, default to help
        if len(argv) == 2:
            argv.append("--help")

        args = parser.parse_args(argv[2:3])
        getattr(dataset_tool, args.command)()

    def learn(self):
        learning_tool = KoverLearningTool()

        parser = argparse.ArgumentParser(usage= \
'''%(prog)s learn <experiment> [<args>]
The most commonly used commands are:
    scm        Learn a conjunction/disjunction model using the Set Covering Machine algorithm.
    tree       Learn a decision tree model using the Classification And Regression Trees algorithm.''')

        parser.add_argument('command', help='The learning experiment to perform',
                            choices=learning_tool.available_commands)

        # If no argument has been specified, default to help
        if len(argv) == 2:
            argv.append("--help")

        args = parser.parse_args(argv[2:3])
        getattr(learning_tool, args.command)()

if __name__ == '__main__':
    CommandLineInterface()
