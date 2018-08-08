#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    Kover: Learn interpretable computational phenotyping models from k-merized genomic data
    Copyright (C) 2015  Alexandre Drouin & Gael Letarte St-Pierre

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
import logging
import numpy as np

from collections import defaultdict
from copy import deepcopy
from functools import partial
from itertools import product
from math import ceil, exp, log as ln, pi, sqrt
from multiprocessing import Pool, cpu_count
from scipy.misc import comb

from ...dataset.ds import KoverDataset
from ..learners.cart import DecisionTreeClassifier, _prune_tree
from ..common.models import CARTModel
from ..common.rules import LazyKmerRuleList, KmerRuleClassifications
from ...utils import _duplicate_last_element, _init_callback_functions, _unpack_binary_bytes_from_ints, _parse_kmer_blacklist
from ..experiments.metrics import _get_binary_metrics, _get_multiclass_metrics


TIEBREAKER_RISK_BLOCK_SIZE = 100000


class BetweenDict(dict):
    """
    The BetweenDict
    By: Joshua Kugler
    Source: http://joshuakugler.com/archives/30-BetweenDict,-a-Python-dict-for-value-ranges.html
    """
    def __init__(self, d = {}):
        for k,v in d.items():
            self[k] = v

    def __getitem__(self, key):
        for k, v in self.items():
            if k[0] <= key < k[1] or (k[0] <= key and k[1] == np.infty) or (k[0] == -np.infty and key < k[1]):
                return v
        raise KeyError("Key '%s' is not between any values in the BetweenDict" % key)

    def __setitem__(self, key, value):
        try:
            if len(key) == 2:
                if key[0] < key[1]:
                    dict.__setitem__(self, (key[0], key[1]), value)
                else:
                    raise RuntimeError('First element of a BetweenDict key '
                                       'must be strictly less than the '
                                       'second element. Got [%.6f, %.6f]' % (key[0], key[1]))
            else:
                raise ValueError('Key of a BetweenDict must be an iterable '
                                 'with length two')
        except TypeError:
            raise TypeError('Key of a BetweenDict must be an iterable '
                             'with length two')

    def __contains__(self, key):
        try:
            return bool(self[key]) or True
        except KeyError:
            return False


def _tiebreaker(best_score_idx, rule_kmer_occurrences):
    """
    The tiebreaker finds the rules in best_score_idx for which the k-mer has the
    largest number of occurrences in the training set. This makes sense because
    it should avoid to select k-mers that isolate a small set of genomes (e.g.,
    a phylogenetic effect) when possible. It should also help find small
    compression sets in the sample compression bound.

    """
    logging.debug("There are %d candidate rules." % best_score_idx.shape[0])
    tie_rule_kmer_occurrences = rule_kmer_occurrences[best_score_idx]
    result = best_score_idx[np.isclose(tie_rule_kmer_occurrences, tie_rule_kmer_occurrences.max())]
    return result


def _split_callback(node, equivalent_rules_idx):
    """
    Called each time a leaf is split into two new leaves. This function is used
    to keep track of the number of equivalent rules at each decision point. It
    simply injects a list of equivalent rule indices into the rule used to split.
    It's kinda hacky, but also really simple.

    """
    logging.debug("Injecting equivalent rules into new tree node")
    node.rule.equivalent_rules_idx = equivalent_rules_idx


def _readdress_tree(tree, rule_new_idx_by_kmer_seq):
    def _readdress(node, kmer_idx):
        if node.rule is not None:
            node.rule.kmer_index = kmer_idx[node.rule.kmer_sequence]
            _readdress(node.left_child, kmer_idx)
            _readdress(node.right_child, kmer_idx)
    new_tree = deepcopy(tree)
    _readdress(new_tree, rule_new_idx_by_kmer_seq)
    return new_tree


def _predictions(decision_tree, kmer_matrix, train_example_idx, test_example_idx, progress_callback=None):
    """
    Makes predictions by loading only the columns of the kmer matrix that are targetted by the model.

    """
    # TODO check if this has the same bug that we fixed for kover SCM
    if progress_callback is None:
        progress_callback = lambda t, p: None
    progress_callback("Testing", 0.0)

    # Case: Standart tree
    if len(decision_tree.rules) > 0:
        model_rules = decision_tree.rules
        kmer_idx_by_rule = np.array([r.kmer_index for r in model_rules])
        kmer_sequence_by_rule = np.array([r.kmer_sequence for r in model_rules])
        sort_by_idx = np.argsort(kmer_idx_by_rule)
        kmer_idx_by_rule = kmer_idx_by_rule[sort_by_idx]
        kmer_sequence_by_rule = kmer_sequence_by_rule[sort_by_idx]
        readdressed_kmer_idx_by_rule = dict((s, i) for i, s in enumerate(kmer_sequence_by_rule))
        readdressed_decision_tree = _readdress_tree(tree=decision_tree, rule_new_idx_by_kmer_seq=readdressed_kmer_idx_by_rule)
        X = np.vstack((_unpack_binary_bytes_from_ints(kmer_matrix[:, idx]) for idx in kmer_idx_by_rule)).T
        train_predictions = readdressed_decision_tree.predict(X[train_example_idx])
        progress_callback("Testing", 1.0 * len(train_example_idx) / (len(train_example_idx) + len(test_example_idx)))
        test_predictions = readdressed_decision_tree.predict(X[test_example_idx])
        progress_callback("Testing", 1.0)

    # Case: the model is just a leaf
    else:
        train_predictions = decision_tree.predict(np.empty((len(train_example_idx), 1)))
        progress_callback("Testing", 1.0 * len(train_example_idx) / (len(train_example_idx) + len(test_example_idx)))
        test_predictions = decision_tree.predict(np.empty((len(test_example_idx), 1)))
        progress_callback("Testing", 1.0)
    return train_predictions, test_predictions


def _bound(train_predictions, train_answers, train_example_idx, model, delta,
           max_genome_size, rule_classifications, n_classes):
    """
    Calculate the value of the sample compression bound for a given decision tree

    """
    logging.debug("Constructing the compression set.")

    # For each rule in the model, find in which genome the corresponding k-mer is present.
    # We build a binary matrix with one row per genome and one column per k-mer
    compression_set = []
    if len(model.rules) > 0:
        # XXX: Otherwise, the model is just one big leaf containing all the training
        #      examples and the compression set is empty.
        presence_by_example = rule_classifications.get_columns([r.kmer_index for r in model.rules])[train_example_idx]

        # Construct the smallest possible compression set (Chvatal greedy approx for minimum set cover)
        # ---
        # While there are still some k-mers to cover
        while presence_by_example.shape[1] != 0:
            # The score is the number of k-mers that a genome contains
            score = presence_by_example.sum(axis=1)
            best_example_relative_idx = np.argmax(score)  # idx of the best training example
            compression_set.append(best_example_relative_idx)

            # Keep only k-mers that were not in the selected genome
            # XXX: no need to remove the selected genome since it doesn't contain any of the remaining k-mers (won't be selected)
            presence_by_example = presence_by_example[:, presence_by_example[best_example_relative_idx] == 0]

    logging.debug("The compression set contains %d examples." % len(compression_set))

    # Compute the bound value
    logging.debug("Computing the bound value.")
    m = float(len(train_answers))  # Number of training examples
    Z_card = float(len(compression_set))  # Number of examples in the compression set
    N_Z = Z_card * max_genome_size  # Number of nucleotides in the compression set
    # Number of errors on examples not in the compression set
    r = float((train_predictions != train_answers).sum() - (train_predictions[compression_set] != train_answers[compression_set]).sum())
    n = float(len(model.rules))

    # Sample-compression bound value
    return 1.0 - exp((-1.0 / (m - Z_card - r)) * (ln(comb(m, Z_card, exact=True)) +
                                                  ln(comb(m - Z_card, r, exact=True)) +
                                                  (n * ln(N_Z) if n > 0 else 0.) +
                                                  (n + 1) * ln(n_classes) +
                                                  ln(comb(2 * n + 1, n, exact=True)) +
                                                  ln(pi**6 *
                                                  (n + 1)**2 *
                                                  (r + 1)**2 *
                                                  (Z_card + 1)**2 /
                                                  (216 * delta))))


def _learn_pruned_tree_bound(hps, dataset_file, split_name, delta, max_genome_size, rule_blacklist):
    """
    Learns a cost-complexity pruned decision tree for a fixed set of hyperparameters and returns an estimate of its
    generalization error.

    The pruning is done based on the bound value.

    Parameters:
    -----------
    dataset_file: str
        The path to the Kover dataset
    split_name: str
        The name of the train/test split to use
    hps: dict
        A dictionnary of hyperparameter values (one value per key)
    rule_blacklist: list
        A list giving the rules to blacklist all the time.

    Returns:
    --------
    hps: dict
        The specified hyperparameter values
    bound_value: float
        The value of the sample-compression bound for the best tree that was grown using the provided HPs
    best_master_tree: float
        A tree that was grown from the entire training set. It is the one with the best generalization performance
        estimated by the sample-compression bound. It has also been pruned to prevent overfitting.

    Notes:
    ------
    The pruning is done using the sample-compression bound for decision trees proposed by Drouin et al. (2017).

    """
    # Open the dataset and load some stuff info into memory
    logging.debug("Loading the kover dataset information")
    dataset = KoverDataset(dataset_file)
    split = dataset.get_split(split_name)
    split.train_genome_idx = split.train_genome_idx[...]
    example_labels = dataset.phenotype.metadata[...]
    n_classes = len(dataset.phenotype.tags)
    rules = LazyKmerRuleList(dataset.kmer_sequences, dataset.kmer_by_matrix_column)
    rule_classifications = KmerRuleClassifications(dataset.kmer_matrix, dataset.genome_count)

    # Initialize the tree to be grown
    master_predictor = DecisionTreeClassifier(criterion=hps["criterion"],
                                              max_depth=hps["max_depth"],
                                              min_samples_split=hps["min_samples_split"],
                                              class_importance=hps["class_importance"])

    # Build an overgrown decision tree on the entire dataset
    logging.debug("Growing the overgrown master tree")
    master_predictor.fit(rules=rules,
                         rule_classifications=rule_classifications,
                         example_idx = {c: split.train_genome_idx[example_labels[split.train_genome_idx] == c]\
                                           for c in range(n_classes)},
                         rule_blacklist=rule_blacklist,
                         tiebreaker=partial(_tiebreaker, rule_kmer_occurrences=rule_classifications.sum_rows(split.train_genome_idx)),
                         level_callback=None,
                         split_callback=_split_callback)

    logging.debug("Pruning the master tree using minimum cost-complexity pruning and the sample-compression bound")
    min_score = np.infty
    min_score_tree = None
    train_answers = example_labels[split.train_genome_idx]
    for alpha, tree in zip(*_prune_tree(master_predictor.decision_tree)):
        train_predictions = _predictions(decision_tree=tree,
                                         kmer_matrix=dataset.kmer_matrix,
                                         train_example_idx=split.train_genome_idx,
                                         test_example_idx=[])[0]

        bound_value = _bound(train_predictions=train_predictions,
                             train_answers=train_answers,
                             train_example_idx=split.train_genome_idx,
                             model=tree,
                             delta=delta,
                             max_genome_size=max_genome_size,
                             rule_classifications=KmerRuleClassifications(dataset.kmer_matrix, dataset.genome_count),
                             n_classes=len(dataset.phenotype.tags))

        if bound_value <= min_score:  # Note: assumes that alphas are sorted in increasing order (so we are always preferring trees that are more pruned)
            min_score = bound_value
            min_score_tree = tree
            hps["pruning_alpha"] = alpha  # Save the best value of alpha
    logging.debug("Pruning completed.")

    # Return the best tree and its error estimate
    return hps, min_score, min_score_tree


def _learn_pruned_tree_cv(hps, dataset_file, split_name, rule_blacklist):
    """
    Learns a cost-complexity pruned decision tree for a fixed set of hyperparameters and returns an estimate of its
    generalization error.

    The pruning is done using k-fold cross-validation.

    Parameters:
    -----------
    dataset_file: str
        The path to the Kover dataset
    split_name: str
        The name of the train/test split to use
    hps: dict
        A dictionnary of hyperparameter values (one value per key)
    rule_blacklist: list
        A dictionnary giving the rules to blacklist all the time.

    Returns:
    --------
    hps: dict
        The specified hyperparameter values
    cv_score: float
        The cross-validation score of the best tree that was grown using the provided HPs
    best_master_tree: float
        A tree that was grown from the entire training set. It is the one with the best generalization performance
        estimated by CV. It has also been pruned to prevent overfitting.

    Notes:
    ------
    The pruning is done using the cross-validation version of the cost-complexity pruning procedure described in:
    Breiman, L., Friedman, J., Stone, C. J., & Olshen, R. A. (1984). Classification and regression trees. CRC press.

    """
    # Open the dataset and load some stuff info into memory
    logging.debug("Loading the kover dataset information")
    dataset = KoverDataset(dataset_file)
    split = dataset.get_split(split_name)
    split.train_genome_idx = split.train_genome_idx[...]
    example_labels = dataset.phenotype.metadata[...]
    n_classes = len(dataset.phenotype.tags)
    rules = LazyKmerRuleList(dataset.kmer_sequences, dataset.kmer_by_matrix_column)
    rule_classifications = KmerRuleClassifications(dataset.kmer_matrix, dataset.genome_count)

    # Initialize the trees to be grown
    logging.debug("Planting seeds")
    fold_predictors = [DecisionTreeClassifier(criterion=hps["criterion"],
                                              max_depth=hps["max_depth"],
                                              min_samples_split=hps["min_samples_split"],
                                              class_importance=hps["class_importance"])
                                              for _ in xrange(len(split.folds))]

    master_predictor = DecisionTreeClassifier(criterion=hps["criterion"],
                                              max_depth=hps["max_depth"],
                                              min_samples_split=hps["min_samples_split"],
                                              class_importance=hps["class_importance"])

    # For each fold, build an overgrown decision tree
    logging.debug("Growing the cross-validation fold trees")
    for i, fold in enumerate(split.folds):
        logging.debug("Growing the tree for fold %d" % (i + 1))
        # Load stuff into memory
        fold.train_genome_idx = fold.train_genome_idx[...]
        fold.test_genome_idx = fold.test_genome_idx[...]

        # Fit the decision tree
        fold_predictors[i].fit(rules=rules,
                               rule_classifications=rule_classifications,
                               example_idx = {c: fold.train_genome_idx[example_labels[fold.train_genome_idx] == c]\
                                                                                            for c in range(n_classes)},
                               rule_blacklist=rule_blacklist,
                               tiebreaker=partial(_tiebreaker, rule_kmer_occurrences=rule_classifications.sum_rows(fold.train_genome_idx)),
                               level_callback=None,
                               split_callback=None)

    # Also build an overgrown decision tree on the entire dataset
    logging.debug("Growing the master tree")
    master_predictor.fit(rules=rules,
                         rule_classifications=rule_classifications,
                         example_idx = {c: split.train_genome_idx[example_labels[split.train_genome_idx] == c]\
                                                                                            for c in range(n_classes)},
                         rule_blacklist=rule_blacklist,
                         tiebreaker=partial(_tiebreaker, rule_kmer_occurrences=rule_classifications.sum_rows(split.train_genome_idx)),
                         level_callback=None,
                         split_callback=_split_callback)

    # Get the pruned master and cross-validation trees
    master_alphas, master_pruned_trees = _prune_tree(master_predictor.decision_tree)
    fold_alphas = []
    fold_pruned_trees = []
    for i in xrange(len(split.folds)):
        alphas, trees = _prune_tree(fold_predictors[i].decision_tree)
        fold_alphas.append(alphas)
        fold_pruned_trees.append(trees)

    # Compute the test risk for all pruned trees of each fold
    fold_scores_by_alpha = []
    for i, fold in enumerate(split.folds):
        fold_test_example_idx = fold.test_genome_idx[...]
        fold_example_labels = example_labels[fold_test_example_idx]
        fold_test_risks = []
        bro = BetweenDict()
        for j, t in enumerate(fold_pruned_trees[i]):
            fold_test_risk = _get_binary_metrics(predictions=_predictions(decision_tree=t,
                                                                          kmer_matrix=dataset.kmer_matrix,
                                                                          train_example_idx=[],
                                                                          test_example_idx=fold_test_example_idx)[1],
                                                 answers=fold_example_labels)["risk"][0]

            fold_test_risks.append(fold_test_risk)
            if j < len(fold_alphas[i]) - 1:
                key = (fold_alphas[i][j], fold_alphas[i][j + 1])
            else:
                key = (fold_alphas[i][j], np.infty)
            bro[key] = fold_test_risk
        fold_scores_by_alpha.append(bro)

    # Prune the master tree based on the CV estimates
    min_score = np.infty
    min_score_tree = None

    logging.debug("The master alphas are: " + str(master_alphas))
    for i, t in enumerate(master_pruned_trees):

        if i < len(master_alphas) - 1:
            geo_mean_alpha_k = sqrt(master_alphas[i] * master_alphas[i + 1])
        else:
            geo_mean_alpha_k = np.infty

        cv_score = np.mean([fold_scores_by_alpha[j][geo_mean_alpha_k] for j in xrange(len(split.folds))])

        if cv_score <= min_score:  # Note: assumes that alphas are sorted in increasing order (so we are always preferring trees that are more pruned)
            min_score = cv_score
            min_score_tree = t
            hps["pruning_alpha"] = geo_mean_alpha_k  # Save the best value of alpha

    # Return the best tree and its error estimate
    return hps, min_score, min_score_tree


def train_tree(dataset_file, split_name, criterion, class_importance, max_depth,
               min_samples_split, rule_blacklist, n_cpu, progress_callback,
               warning_callback, error_callback, hp_search_func, hp_search_type):
    """
    Train a decision tree classifier with the best hyperparameter values, which
    are selected according to hp_search_func.

    """
    # Find the best combination of hyperparameters
    n_hp_combinations = len(criterion) * len(class_importance) * len(max_depth) * len(min_samples_split)
    logging.debug("There are %d hyperparameter combinations to try." % n_hp_combinations)

    logging.debug("Using %d CPUs." % n_cpu)
    pool = Pool(n_cpu)
    _hp_eval_func = partial(hp_search_func, dataset_file=dataset_file, split_name=split_name, rule_blacklist=rule_blacklist)
    best_hps = None
    best_score = np.infty
    best_master_tree = None

    n_completed = 0.0
    progress_callback(hp_search_type.title(), 0.0)
    for hps, score, master_tree in pool.imap_unordered(_hp_eval_func,
                                                ({"criterion": hps[0],
                                                  "class_importance": hps[1],
                                                  "max_depth": hps[2],
                                                  "min_samples_split": hps[3]}
                                                  for hps in product(criterion, class_importance, max_depth,
                                                                     min_samples_split))):

        # XXX: master_tree is the pruned decision tree learned on the entire training set (no need to retrain later)
        n_completed += 1
        progress_callback(hp_search_type.title(), n_completed / n_hp_combinations)
        if score < best_score:
            best_hps = hps
            best_score = score
            best_master_tree = master_tree
        elif np.isclose(score, best_score):
            master_tree_length = len(master_tree)
            best_master_tree_length = len(best_master_tree)
            # XXX: In case of a tie, use the following rules (in order of precedence):
            #      1. Pick the smallest tree
            #      2. Pick the tree with the least variance in the class importances
            #      Note: max_depth and min_samples_split are not specified as they are implied by rule 1
            if (master_tree_length < best_master_tree_length) or \
               (master_tree_length == best_master_tree_length and np.var(hps["class_importance"].values()) < np.var(best_hps["class_importance"].values())):
                best_hps = hps
                best_master_tree = best_master_tree
                best_score = score


    return best_score, best_hps, best_master_tree
    

def _find_rule_blacklist(dataset_file, kmer_blacklist_file, warning_callback):
    """
    Finds the index of the rules that must be blacklisted.
    """
    dataset = KoverDataset(dataset_file)

    # Find all rules to blacklist
    rule_blacklist = []
    if kmer_blacklist_file is not None:
        kmers_to_blacklist = _parse_kmer_blacklist(kmer_blacklist_file, dataset.kmer_length)

        if kmers_to_blacklist:
            # XXX: the k-mers are assumed to be upper-cased in the dataset
            kmer_sequences = dataset.kmer_sequences[...].tolist()
            kmer_by_matrix_column = dataset.kmer_by_matrix_column[...].tolist() # XXX: each k-mer is there only once (see wiki)
            n_kmers = len(kmer_sequences)

            kmers_not_found = []
            for k in kmers_to_blacklist:
                k = k.upper()
                try:
                    rule_blacklist.append(kmer_by_matrix_column.index(kmer_sequences.index(k))) # XXX: We only consider presence rules
                except ValueError:
                    kmers_not_found.append(k)

            if(len(kmers_not_found) > 0):
                warning_callback("The following kmers could not be found in the dataset: " + ", ".join(kmers_not_found))

    return rule_blacklist


def learn_CART(dataset_file, split_name, criterion, max_depth, min_samples_split,
               class_importance, bound_delta, bound_max_genome_size, kmer_blacklist_file,
               parameter_selection, n_cpu, authorized_rules,
               progress_callback=None, warning_callback=None, error_callback=None):
    """
    Cross-validate the best hyper-parameters (criterion, max_depth, min_samples_split and class_importance)
    to grow a pruned decision tree.

    """
    # Initialize callback functions
    warning_callback, error_callback, progress_callback = _init_callback_functions(warning_callback, error_callback,
                                                                                   progress_callback)
    logging.debug("Searching for blacklisted rules.")
    rule_blacklist = _find_rule_blacklist(dataset_file=dataset_file,
                                          kmer_blacklist_file=kmer_blacklist_file,
                                          warning_callback=warning_callback)
                                          
    # Load the dataset info
    dataset = KoverDataset(dataset_file)

    # Check and initialize (hyper)parameters
    if n_cpu is None:
        n_cpu = cpu_count()
    criterion = np.unique(criterion)
    class_importance = np.unique(class_importance)
    max_depth = np.unique(max_depth)
    min_samples_split = np.unique(min_samples_split)

    if parameter_selection == "bound":
        func = partial(_learn_pruned_tree_bound, delta=bound_delta,
                                                 max_genome_size=bound_max_genome_size)
        best_hp_score, best_hps, best_master_tree = \
            train_tree(hp_search_func=func,
                       hp_search_type="bound selection",
                       dataset_file=dataset_file,
                       split_name=split_name,
                       criterion=criterion,
                       class_importance=class_importance,
                       max_depth=max_depth,
                       min_samples_split=min_samples_split,
                       rule_blacklist=rule_blacklist,
                       n_cpu=n_cpu,
                       progress_callback=progress_callback,
                       warning_callback=warning_callback,
                       error_callback=error_callback)

    elif parameter_selection == "cv":
        n_folds = len(dataset.get_split(split_name).folds)
        if n_folds < 1:
            error_callback(Exception("Cross-validation cannot be performed on a split with no folds."))

        best_hp_score, best_hps, best_master_tree = \
            train_tree(hp_search_func=_learn_pruned_tree_cv,
                       hp_search_type="cross-validation",
                       dataset_file=dataset_file,
                       split_name=split_name,
                       criterion=criterion,
                       class_importance=class_importance,
                       max_depth=max_depth,
                       min_samples_split=min_samples_split,
                       rule_blacklist=rule_blacklist,
                       n_cpu=n_cpu,
                       progress_callback=progress_callback,
                       warning_callback=warning_callback,
                       error_callback=error_callback)

    else:
        error_callback(ValueError("Unknown hyperparameter selection strategy specified."))

    # Open the dataset and load some split info into memory
    logging.debug("Opening the Kover dataset and loading split information into memory")
    dataset = KoverDataset(dataset_file)
    split = dataset.get_split(split_name)
    split.train_genome_idx = split.train_genome_idx[...]
    split.test_genome_idx = split.test_genome_idx[...]
    example_labels = dataset.phenotype.metadata[...]
    phenotype_tags = dataset.phenotype.tags[...]

    # Using the best hyperparameters, compute predictions and metrics
    train_predictions, test_predictions = _predictions(decision_tree=best_master_tree,
                                                       kmer_matrix=dataset.kmer_matrix,
                                                       train_example_idx=split.train_genome_idx,
                                                       test_example_idx=split.test_genome_idx,
                                                       progress_callback=progress_callback)

    train_answers = example_labels[split.train_genome_idx]
    test_answers = example_labels[split.test_genome_idx]

    if dataset.classification_type == "binary":
        train_metrics = _get_binary_metrics(train_predictions, train_answers)
    else:
        train_metrics = _get_multiclass_metrics(train_predictions, train_answers, len(phenotype_tags))
    if len(split.test_genome_idx) > 0:
        if dataset.classification_type == "binary":
            test_metrics = _get_binary_metrics(test_predictions, test_answers)
        else:
            test_metrics = _get_multiclass_metrics(test_predictions, test_answers, len(phenotype_tags))
    else:
        test_metrics = None

    # Get the idx of the training/testing examples that are correctly/incorrectly classified by the model
    classifications = defaultdict(list)
    classifications["train_correct"] = dataset.genome_identifiers[split.train_genome_idx[train_predictions == \
                                                train_answers].tolist()].tolist() if train_metrics["risk"][0] < 1.0 else []
    classifications["train_errors"] = dataset.genome_identifiers[split.train_genome_idx[train_predictions != \
                                                train_answers].tolist()].tolist() if train_metrics["risk"][0] > 0 else []
    if len(split.test_genome_idx) > 0:
        classifications["test_correct"] = dataset.genome_identifiers[split.test_genome_idx[test_predictions == \
                                                test_answers].tolist()].tolist() if test_metrics["risk"][0] < 1.0 else []
        classifications["test_errors"] = dataset.genome_identifiers[split.test_genome_idx[test_predictions != \
                                                test_answers].tolist()].tolist() if test_metrics["risk"][0] > 0 else []

    best_model = CARTModel(class_tags=phenotype_tags)
    best_model.decision_tree = best_master_tree

    # Extract all the equivalent rules for the nodes in the model
    rules = LazyKmerRuleList(dataset.kmer_sequences, dataset.kmer_by_matrix_column)
    model_equivalent_rules = {r: [rules[i] for i in r.equivalent_rules_idx] for r in best_master_tree.rules}

    # Extract the importance of each node in the model and normalize it
    rule_importance_sum = float(sum(r.importance for r in best_master_tree.rules))
    rule_importances = {r: r.importance / rule_importance_sum for r in best_master_tree.rules}

    return best_hps, best_hp_score, train_metrics, test_metrics, best_model,\
           rule_importances, model_equivalent_rules, \
           classifications
