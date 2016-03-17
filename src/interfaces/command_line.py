#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
	Kover: Learn interpretable computational phenotyping models from k-merized genomic data
	Copyright (C) 2015  Alexandre Drouin

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

from sys import argv

KOVER_DESCRIPTION = "Kover: Learn interpretable computational phenotyping models from k-merized genomic data"
VERSION = "0.1"


class KoverDatasetTool(object):
    def __init__(self):
        self.available_commands = ['create', 'info', 'split']

    def create(self):
        parser = argparse.ArgumentParser(prog="kover dataset create",
                                         description='Creates a Kover dataset from genomic data and optionally phenotypic metadata')
        parser.add_argument('--genome-type', choices=['tsv'], help='The format in which the genomic data is '
                            'provided. See documentation for details.', required=True)
        parser.add_argument('--genome-source', help='The genomic data.',
                            required=True)
        parser.add_argument('--phenotype-name', help='An informative name that is assigned to the phenotypic metadata.')
        parser.add_argument('--phenotype-metadata', help='A file containing the phenotypic metadata.')
        parser.add_argument('--output', help='The Kover dataset to be created.', required=True)
        parser.add_argument('--compression', type=int, help='The gzip compression level (0 - 9). 0 means no compression'
                            '. The default value is 4.', default=4)
        parser.add_argument('-x', '--progress', help='Shows a progress bar for the execution.', action='store_true')
        parser.add_argument('-v', '--verbose', help='Sets the verbosity level.', default=False, action='store_true')

        # If no argument has been specified, default to help
        if len(argv) == 3:
            argv.append("--help")

        args = parser.parse_args(argv[3:])

        # Input validation logic
        if (args.phenotype_name is not None and args.phenotype_metadata is None) or (
                        args.phenotype_name is None and args.phenotype_metadata is not None):
            print "Error: The phenotype name and metadata file must be specified."
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

        if args.genome_type == "tsv":
            from kover.core.dataset.create import from_tsv

            from_tsv(tsv_path=args.genome_source,
                     output_path=args.output,
                     phenotype_name=args.phenotype_name,
                     phenotype_metadata_path=args.phenotype_metadata,
                     gzip=args.compression,
                     progress_callback=progress)

        if args.progress:
            progress_vars["pbar"].finish()

    def info(self):
        # TODO: list the available splits
        parser = argparse.ArgumentParser(prog="kover dataset info",
                                         description='Prints information about the content of a dataset')
        parser.add_argument('--dataset', help='The Kover dataset for which you require information.', required=True)
        parser.add_argument('--all', help='Prints all the available information.', action='store_true')
        parser.add_argument('--genome-type', help='Prints the type of genomic data that was used to create the dataset.',
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
        parser.add_argument('--phenotype-name', help='Prints the identifier that was assigned to the phenotype.',
                            action='store_true')
        parser.add_argument('--phenotype-metadata', help='Prints the path of the file from which the phenotypic metadata'
                            ' was acquired.', action='store_true')
        parser.add_argument('--splits', help='Prints the lists of splits of the dataset that are available for learning.',
                            action='store_true')
        parser.add_argument('--uuid', help='Prints the unique identifier of the dataset.', action='store_true')
        parser.add_argument('--compression', help='Print the data compression options of the dataset.',
                            action='store_true')

        # If no argument has been specified, default to help
        if len(argv) == 3:
            argv.append("--help")

        args = parser.parse_args(argv[3:])

        # Package imports
        from kover.core.dataset import KoverDataset

        dataset = args.dataset
        dataset = KoverDataset(dataset)
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
        if args.phenotype_name or args.all:
            print "Phenotype name:", dataset.phenotype.name
            print
        if args.phenotype_metadata or args.all:
            if dataset.phenotype.name != "NA":
                print "Phenotype metadata source:", dataset.phenotype.metadata_source
            else:
                print "No phenotype metadata."
            print
        if args.compression or args.all:
            print "Compression:", dataset.compression
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
                            'training the learning algorithm (default is 0.5). Alternatively, you can specify which '
                            'genomes to use for training and testing by using --train-ids and --test-ids.', default=0.5)
        parser.add_argument('--train-ids', type=str, nargs='+', help='The identifiers of the genomes used to train the '
                            'learning algorithm. If you provide a value for this argument, you must also provide a '
                             'value for --test-ids.')
        parser.add_argument('--test-ids', type=str, nargs='+', help='The identifiers of the genomes used to evaluate '
                            'the accuracy of the model generated. If you provide a value for this argument, you must '
                            'also provide a value for --train-ids.')
        parser.add_argument('--folds', type=int,
                            help='The number of k-fold cross-validation folds to create (default is 0 for none, '
                                 'the minimum value is 2). Folds are required for using k-fold cross-validation '
                                 'in \'kover learn\'.', default=0)
        parser.add_argument('--random-seed', type=int, help='A random seed used for randomly splitting the data. A '
                            'specific seed will always lead to the same split. If not provided, it is set randomly.')
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
        from kover.core.dataset.split import split_with_ids, split_with_proportion
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
                           train_ids=args.train_ids,
                           test_ids=args.test_ids,
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
                    Copyright (C) 2015  Alexandre Drouin

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
                print "kover version %s" % VERSION
            elif args.cite:
                print \
                    """If you used Kover in your work, please cite our paper:

Drouin, A., Giguère, S., Déraspe, M., Marchand, M., Tyers, M., Loo, V. G., Bourgault, A. M., Laviolette, F. & Corbeil, J. (2015). Predictive computational phenotyping and biomarker discovery using reference-free genome comparisons. submitted.
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
        parser = argparse.ArgumentParser(prog='kover learn', description='Learn a model from data')
        parser.add_argument('--dataset', help='The Kover dataset from which to learn.', required=True)
        parser.add_argument('--split', help='The identifier of the split of the dataset to use for learning.',
                            required=True)
        parser.add_argument('--model-type', choices=['conjunction', 'disjunction'], nargs='+',
                            help='Hyperparameter: The type of model to learn, conjunction (logical-AND) or disjunction '
                                 '(logical-OR). You can specify multiple space separated values.',
                            required=True)
        parser.add_argument('--p', type=float, nargs='+',
                            help='Hyperparameter: The value of the trade-off parameter in the rule scoring criterion. '
                                 'You can specify multiple space separated values.',
                            required=True)
        parser.add_argument('--max-rules', type=int, help='The maximum number of rules that can be included in the '
                            'model.', required=True)
        parser.add_argument('--max-equiv-rules', type=int, help='The maximum number of equivalent rules to report for '
                            'each rule in the model. This only affects model interpretation. Use the default unless you'
                            ' expect that the rules in the model will be equivalent to more than 10000 other rules.',
                            default=10000)
        parser.add_argument('--hp-choice', choices=['bound', 'cv', 'none'],
                            help='The strategy used to select the best values for hyperparameters. The default is '
                            'k-fold cross-validation, where k is the number of folds defined in the split. Other '
                            'strategies, such as bound selection are available. Using none selects the first value '
                            'specified for each hyperparameter.', default='cv')
        parser.add_argument('--bound-delta', type=float, help='The probabilistic bound on the error rate will be valid '
                            'with probability 1-delta. The default value is 0.05.', default=0.05)
        parser.add_argument('--bound-max-genome-size', type=int, help='The maximum size, in base pairs, of any '
                            'genome in the dataset. If you are unsure about this value, you should use an '
                            'overestimation. This will only affect the tightness of the bound on the error rate. By '
                            'default number of k-mers in the dataset is used.')
        parser.add_argument('--random-seed', type=int, help='The random seed used for any random operation. '
                            'Set this if only if you require that the same random choices are made between repeats.')
        parser.add_argument('--n-cpu', type=int, help='The number of CPUs used to select the hyperparameter values. '
                            'Make sure your computer has enough RAM to handle multiple simultaneous trainings of the '
                            'algorithm and that your storage device will not be a bottleneck (simultaneous reading).',
                            default=1)
        parser.add_argument('--output-dir', help='The directory in which to store Kover\'s output. It will be created if '
                            'it does not exist.', default='.')
        parser.add_argument('-x', '--progress', help='Shows a progress bar for the execution.', action='store_true')
        parser.add_argument('-v', '--verbose', help='Sets the verbosity level.', default=False, action='store_true')

        # If no argument has been specified, default to help
        if len(argv) == 2:
            argv.append("--help")

        args = parser.parse_args(argv[2:])

        # Package imports
        from datetime import timedelta
        from json import dump as json_dump
        from kover.core.dataset import KoverDataset
        from kover.core.learning.experiment import learn
        from os import mkdir
        from os.path import abspath, exists, join
        from progressbar import Bar, Percentage, ProgressBar, Timer
        from time import time

        # Input validation
        dataset = KoverDataset(args.dataset)
        dataset_kmer_count = dataset.kmer_count
        # - Check that the split exists
        try:
            dataset.get_split(args.split)
        except:
            print "Error: The split (%s) does not exist in the dataset. Use 'kover dataset split' to create it." % args.split
            exit()
        # - Must have at least 2 folds to perform cross-validation
        if args.hp_choice == "cv" and len(dataset.get_split(args.split).folds) < 2:
            print "Error: The split must contain at least 2 folds in order to perform cross-validation. " \
                  "Use 'kover dataset split' to create folds."
            exit()
        del dataset

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

        start_time = time()
        best_hp, best_hp_score, \
        train_metrics, test_metrics, \
        model, rule_importances, \
        equivalent_rules, \
        classifications = learn(dataset_file=args.dataset,
                                split_name=args.split,
                                model_type=args.model_type,
                                p=args.p,
                                max_rules=args.max_rules,
                                max_equiv_rules=args.max_equiv_rules,
                                bound_delta = args.bound_delta,
                                bound_max_genome_size = args.bound_max_genome_size
                                                        if args.bound_max_genome_size is not None
                                                        else dataset_kmer_count,
                                parameter_selection=args.hp_choice,
                                n_cpu=args.n_cpu,
                                random_seed=args.random_seed,
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
        report += "Phenotype: %s\n" % dataset.phenotype.name.title()
        report += "Split: %s\n" % args.split
        report += "Number of genomes used for training: %d (Group 1: %d, Group 0: %d)\n" % (
            len(split.train_genome_idx),
            (dataset.phenotype.metadata[split.train_genome_idx] == 1).sum(),
            (dataset.phenotype.metadata[split.train_genome_idx] == 0).sum())
        report += "Number of genomes used for testing: %d (Group 1: %d, Group 0: %d)\n" % (
            len(split.test_genome_idx),
            (dataset.phenotype.metadata[split.test_genome_idx] == 1).sum() if len(split.test_genome_idx) > 0 else 0,
            (dataset.phenotype.metadata[split.test_genome_idx] == 0).sum() if len(split.test_genome_idx) > 0 else 0)
        report += "Number of k-mers: %d\n" % dataset.kmer_count
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
        # Print the training set metrics
        report += "Metrics (training data)\n" + "-" * 23 + "\n"
        for key, alias in metric_aliases:
            report += "%s: %s\n" % (str(alias), str(round(train_metrics[key][0], 5)))
        report += "\n"
        # Print the testing set metrics
        if test_metrics is not None:
            report += "Metrics (testing data)\n" + "-" * 22 + "\n"
            for key, alias in metric_aliases:
                report += "%s: %s\n" % (str(alias), str(round(test_metrics[key][0], 5)))
            report += "\n"
        # Print the bound
        report += "Metrics (Probabilistic bound)\n" + "-" * 22 + "\n"
        report += "With probability %.1f%%, the model will have a maximum error rate of %.1f%% on any unseen genome." % ((1.0 - args.bound_delta) * 100, train_metrics["bound"] * 100)
        report += "\n"
        report += "\n"

        report += "Model (%s - %d rules):\n" % (model.type.title(), len(model)) + "-" * (
            18 + len(model.type) + len(str(len(model)))) + "\n"
        report += ("\n%s\n" % ("AND" if model.type == "conjunction" else "OR")).join(
            ["%s [Importance: %.2f, %d equivalent rules]" % (str(rule), importance, len(equivalent_rules[i])) for
             i, (rule, importance) in enumerate(zip(model, rule_importances))])
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

        # Save model (also equivalent rules) [Fasta]
        with open(join(args.output_dir, 'model.fasta'), "w") as f:
            for i, (rule, importance) in enumerate(zip(model, rule_importances)):
                f.write(">rule-%d %s, importance: %.2f\n%s\n\n" % (i + 1, rule.type, importance, rule.kmer_sequence))
                with open(join(args.output_dir, "model_rule_%i_equiv.fasta" % (i + 1)), "w") as f_equiv:
                    f_equiv.write("\n\n".join(
                        [">rule-%d-equiv-%d,%s\n%s" % (i + 1, j + 1, rule.type, rule.kmer_sequence)
                         for j, rule in enumerate(equivalent_rules[i])]))

if __name__ == '__main__':
    CommandLineInterface()
