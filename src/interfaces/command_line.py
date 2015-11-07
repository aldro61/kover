#!/usr/bin/env python

import argparse
import logging

from sys import argv

class KoverDatasetTool(object):
    def __init__(self):
        self.available_commands = ['create', 'info', 'split']

    def create(self):
        parser = argparse.ArgumentParser(prog="kover dataset create", description='Creates a Kover dataset from genomic data and optionally phenotypic metadata.')
        parser.add_argument('-t', '--genome-type', choices=['tsv'], help='The type of source for the genomic data.', required=True)
        parser.add_argument('-g', '--genome-source', help='The source of the genomic data.', required=True)
        parser.add_argument('-o', '--output', help='The output Kover dataset.', required=True)
        parser.add_argument('-p', '--phenotype-name', help='The name of the phenotype.', required=False)
        parser.add_argument('-m', '--phenotype-metadata', help='A file containing the metadata.', required=False)
        parser.add_argument('-z', '--compression', type=int, help='The gzip compression level (0 - 9). 0 means no compression. The default value is 4.', default=4, required=False)
        parser.add_argument('-x', '--progress', help='Shows a progress bar for the execution.', action='store_true', required=False)
        parser.add_argument('-v', '--verbose', help='Sets the verbosity level.', default=False, action='store_true', required=False)
        args = parser.parse_args(argv[3:])

        # Input validation logic
        if (args.phenotype_name is not None and args.phenotype_metadata is None) or (args.phenotype_name is None and args.phenotype_metadata is not None):
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
                    progress_vars["pbar"] = ProgressBar(widgets=['%s: ' % task_name, Percentage(), Bar(), Timer()], maxval=1.0)
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
        #TODO: list the available splits
        parser = argparse.ArgumentParser(prog="kover dataset info", description='Prints information about a dataset.')
        parser.add_argument('-d', '--dataset', help='A Kover dataset.', required=True)
        parser.add_argument('-a', '--all', help='Print all the available information.', action='store_true', required=False)
        parser.add_argument('-t', '--genome-type', help='Print the type of source from which Kover obtained the genomic data', action='store_true', required=False)
        parser.add_argument('-g', '--genome-source', help='Print the source from which Kover obtained the genomic data.', action='store_true', required=False)
        parser.add_argument('-i', '--genome-ids', help='Print the genome identifers.', action='store_true', required=False)
        parser.add_argument('-n', '--genome-count', help='Print the number of genomes.', action='store_true', required=False)
        parser.add_argument('-k', '--kmers', help='Print the k-mer sequences in fasta format.', action='store_true', required=False)
        parser.add_argument('-l', '--kmer-len', help='Print the k-mer length.', action='store_true', required=False)
        parser.add_argument('-c', '--kmer-count', help='Print the k-mer count.', action='store_true', required=False)
        parser.add_argument('-p', '--phenotype-name', help='Print the phenotype\'s name.', action='store_true', required=False)
        parser.add_argument('-m', '--phenotype-metadata', help='Print the file from which the metadata was extracted.', action='store_true', required=False)
        parser.add_argument('-u', '--uuid', help='Print the unique identifier of the Kover dataset (data independent).')
        parser.add_argument('-z', '--compression', help='Print the compression options.', action='store_true', required=False)
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
                print ">k%d" % (i+1)
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
            print "Phenotype metadata source:", dataset.phenotype.metadata_source
            print
        if args.compression or args.all:
            print "Compression:", dataset.compression
            print

    def split(self):
        parser = argparse.ArgumentParser(prog="kover dataset split", description='Splits a kover dataset file into a training set, a testing set and optionally cross-validation folds.')
        parser.add_argument('-d', '--dataset', help='The Kover dataset to split.', required=True)
        parser.add_argument('-i', '--id', help='The identifier of the split.', required=True)
        parser.add_argument('-t', '--train-size', type=float, help='The proportion of the data used for training (default is 0.5).', default=0.5, required=False)
        parser.add_argument('-r', '--random-seed', type=int, help='The random seed used for the split (If not provided a random value between 0 and 4294967295 will be used).', required=False)
        parser.add_argument('-f', '--folds', type=int, help='The number of cross-validation folds (default is 0 for none, the minimum value is 2).', default=0, required=False)
        parser.add_argument('-v', '--verbose', help='Sets the verbosity level.', default=False, action='store_true', required=False)
        parser.add_argument('-x', '--progress', help='Shows a progress bar for the execution.', action='store_true', required=False)
        args = parser.parse_args(argv[3:])

        # Input validation
        if args.folds == 1:
            print "Error: The number of cross-validation folds must be 0 or >= 2."
            exit()

        # Package imports
        from kover.core.dataset.split import split
        from progressbar import Bar, Percentage, ProgressBar, Timer
        from random import randint

        if args.verbose:
            logging.basicConfig(level=logging.DEBUG,
                                format="%(asctime)s.%(msecs)d %(levelname)s %(module)s - %(funcName)s: %(message)s")

        if args.randomseed is None:
            args.randomseed = randint(0, 4294967295)

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

        split(input=args.dataset,
              identifier=args.id,
              train_prop=args.trainsize,
              random_seed=args.randomseed,
              n_folds=args.folds,
              progress_callback=progress)

        if args.progress:
            progress_vars["pbar"].finish()


class CommandLineInterface(object):

    def __init__(self):
        self.available_commands = ['dataset', 'learn']

        parser = argparse.ArgumentParser(
            #description='Kover -- Learn from kmers',
            usage='''kover <command> [<args>]

The most commonly used commands are:
    dataset     Dataset manipulation tools
''')
        parser.add_argument('command', help='Command to run')
        # parse_args defaults to [1:] for args, but you need to
        # exclude the rest of the args too, or validation will fail
        args = parser.parse_args(argv[1:2])
        if not hasattr(self, args.command):
            print 'kover: \'%s\' is not a kover command. See \'%s --help\'.' % (args.command, argv[0])

            parser.print_help()
            exit(1)
        # use dispatch pattern to invoke method with same name
        getattr(self, args.command)()

    def dataset(self):
        dataset_tool = KoverDatasetTool()

        parser = argparse.ArgumentParser(description='Dataset manipulation utilities',
                                         usage='''%(prog)s dataset <command> [<args>]
The most commonly used commands are:
    create     Creates a Kover datasets from genomic data.
    split      Splits a kover dataset file into a training set, a testing set and optionally cross-validation folds.
    info       Prints information about a dataset''')

        parser.add_argument('command', help='The dataset manipulation to perform',
                            choices=dataset_tool.available_commands)
        args = parser.parse_args(argv[2:3])
        getattr(dataset_tool, args.command)()

    def learn(self):
        parser = argparse.ArgumentParser(prog='kover learn', description='Learn a model from data')
        parser.add_argument('-d', '--dataset', help='The Kover dataset to learn from.', required=True)
        parser.add_argument('-s', '--split', help='The identifier of the split of the dataset that must be learnt from.', required=True)
        parser.add_argument('-t', '--model-type', choices=['conjunction', 'disjunction'], nargs='+', help='Hyperparameter: The type of model (conjunction or disjunction) to learn. Single value or multiple space separated values.', required=True)
        parser.add_argument('-p', '--p', type=float, nargs='+', help='Hyperparameter: The value of the trade-off used to score the rules. Single value or multiple space separated values.', required=True)
        parser.add_argument('-m', '--max-rules', type=int, help='The maximum number of rules to include in a model.', required=True)
        parser.add_argument('-c', '--hp-choice', choices=['bound', 'cv', 'none'], help='The strategy used to select the hyperparameter values.', default='cv', required=False)
        parser.add_argument('-P', '--n-cpu', type=int, help='The number of CPUs used to select the hyperparameter values. Use 1 unless you are using a parallel file system (see documentation).', default=1, required=False)
        parser.add_argument('-o', '--output-dir', help='The directory in which to store Kover\'s output. Will be created if it does not exist.', default='.', required=False)
        parser.add_argument('-x', '--progress', help='Shows a progress bar for the execution.', action='store_true', required=False)
        parser.add_argument('-v', '--verbose', help='Sets the verbosity level.', default=False, action='store_true', required=False)
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
                    progress_vars["pbar"] = ProgressBar(widgets=['%s: ' % task_name, Percentage(), Bar(), Timer()], maxval=1.0)
                    progress_vars["pbar"].start()
                else:
                    progress_vars["pbar"].update(p)
        else:
            progress = None

        start_time = time()
        best_hp, best_hp_score, \
        train_metrics, test_metrics, model = learn(dataset_file=args.dataset,
                                                   split_name=args.split,
                                                   model_type=args.model_type,
                                                   p=args.p,
                                                   max_rules=args.max_rules,
                                                   parameter_selection=args.hp_choice,
                                                   n_cpu=args.n_cpu,
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
        report += "Phenotype: %s\n" % dataset.phenotype.name.title()
        report += "Split: %s\n" % args.split
        report += "Number of genomes used for training: %d (Positive: %d, Negative: %d)\n" % (len(dataset.get_split(args.split).train_genome_idx), (dataset.phenotype.metadata[dataset.get_split(args.split).train_genome_idx] == 1).sum(), (dataset.phenotype.metadata[dataset.get_split(args.split).train_genome_idx] == 0).sum())
        report += "Number of genomes used for testing: %d (Positive: %d, Negative: %d)\n" % (len(dataset.get_split(args.split).test_genome_idx), (dataset.phenotype.metadata[dataset.get_split(args.split).test_genome_idx] == 1).sum(), (dataset.phenotype.metadata[dataset.get_split(args.split).test_genome_idx] == 0).sum())
        report += "\n"
        report += "Hyperparameter Values:\n" + "-" * 22 + "\n"
        if args.hp_choice == "cv":
            report += "Selection strategy: %d-fold cross-validation (score = %.5f)\n" % (len(dataset.get_split(args.split).folds), best_hp_score)
        elif args.hp_choice == "bound":
            report += "Selection strategy: bound selection (score = %.5f)\n" % best_hp_score
        else:
            report += "Selection strategy: No selection\n"
        report += "Model type: %s\n" % best_hp["model_type"]
        report += "p: %f\n" % best_hp["p"]
        report += "Maximum number of rules: %d\n" % best_hp["max_rules"]
        report += "\n"
        report += "Metrics\n" + "-"*7 + "\n"
        for key, alias in metric_aliases:
            report += "%s: %s\n" % (str(alias), str(round(test_metrics[key][0], 5)))
        report += "\n"
        report += "Model (%s - %d rules):\n" % (model.type.title(), len(model)) + "-"*(18 + len(model.type) + len(str(len(model)))) + "\n"
        report += (("\n%s\n" % ("AND" if model.type == "conjunction" else "OR"))).join([str(rule) for rule in model])
        report += "\n"

        print report

        if not exists(args.output_dir):
            mkdir(args.output_dir)

        # Write report
        with open(join(args.output_dir, "report.txt"), "w") as f:
            f.write(report)

        # Write metrics (not user friendly) [json]
        results = {"data": {"uuid": dataset.uuid,
                            "split": args.split},
                   "cv": {"best_hp": best_hp,
                          "best_hp_score": best_hp_score,
                          "strategy": args.hp_choice},
                   "metrics": {
                       "train": train_metrics,
                       "test": test_metrics},
                   "running_time": running_time.seconds}

        with open(join(args.output_dir, 'results.json'), 'w') as f:
            json_dump(results, f)

        # Save model (also equivalent rules) [json]
        with open(join(args.output_dir, 'model.fasta'), "w") as f:
            for i, rule in enumerate(model):
                f.write(">rule-%d %s\n%s\n\n" % (i + 1, rule.type, rule.kmer_sequence))
        #TODO: save equivalent rules

if __name__ == '__main__':
    CommandLineInterface()