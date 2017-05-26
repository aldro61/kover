#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
	Kover: Learn interpretable computational phenotyping models from k-merized genomic data
	Copyright (C) 2015  Alexandre Drouin & Gaël Letarte St-Pierre

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

from ..dataset.ds import KoverDataset
from .learners.cart import DecisionTreeClassifier, _prune_tree
from .common.rules import LazyKmerRuleList, KmerRuleClassifications
from ..utils import _duplicate_last_element, _init_callback_functions, _unpack_binary_bytes_from_ints

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
            if k[0] <= key < k[1]:
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
            
def _get_metrics(predictions, answers):
    if len(predictions.shape) == 1:
        predictions = predictions.reshape(1, -1)
    metrics = defaultdict(list)
    for i in xrange(predictions.shape[0]):
        p = predictions[i]
        risk = 1.0 * len(p[p != answers]) / len(answers)
        tp = len(np.where(p[answers == 1] == 1)[0])
        fp = len(np.where(p[answers == 0] == 1)[0])
        tn = len(np.where(p[answers == 0] == 0)[0])
        fn = len(np.where(p[answers == 1] == 0)[0])
        precision = 1.0 * tp / (tp + fp) if (tp + fp) != 0 else -np.infty
        sensitivity = recall = 1.0 * tp / (tp + fn) if (tp + fn) != 0 else -np.infty
        specificity = 1.0 * tn / (fp + tn) if (fp + tn) != 0 else -np.infty
        f1_score = 2.0 * precision * recall / (precision + recall) if (precision + recall) > 0.0 else -np.infty
        metrics["risk"].append(risk)
        metrics["tp"].append(tp)
        metrics["fp"].append(fp)
        metrics["tn"].append(tn)
        metrics["fn"].append(fn)
        metrics["precision"].append(precision)
        metrics["sensitivity"].append(sensitivity)
        metrics["recall"].append(recall)
        metrics["specificity"].append(specificity)
        metrics["f1_score"].append(f1_score)
    return metrics
    
def _tiebreaker(best_score_idx, rule_risks):
    """
    The tiebreaker finds the rules in best_score_idx that have the smallest empirical risk.

    """
    logging.debug("There are %d candidate rules." % best_score_idx.shape[0])
    tie_rule_risks = rule_risks[best_score_idx]
    result = best_score_idx[tie_rule_risks == tie_rule_risks.min()]
    return result
    
def _readdress_model(model, rule_new_idx_by_kmer_seq):
    def _readdress(node, kmer_idx):
        if node.rule is not None:
            node.rule.kmer_index = kmer_idx[node.rule.kmer_sequence]
            _readdress(node.left_child, kmer_idx)
            _readdress(node.right_child, kmer_idx)
    new_model = deepcopy(model)
    _readdress(new_model, rule_new_idx_by_kmer_seq)
    return new_model
    
def _predictions(model, kmer_matrix, train_example_idx, test_example_idx, progress_callback=None):
    """
    Makes predictions by loading only the columns of the kmer matrix that are targetted by the model.
    """
    # TODO check if this has the same bug that we fixed for kover SCM
    if progress_callback is None:
        progress_callback = lambda t, p: None
    progress_callback("Testing", 0.0)
    if len(model.rules) > 0:
        model_rules = model.rules
        kmer_idx_by_rule = np.array([r.kmer_index for r in model_rules])
        kmer_sequence_by_rule = np.array([r.kmer_sequence for r in model_rules])
        sort_by_idx = np.argsort(kmer_idx_by_rule)
        kmer_idx_by_rule = kmer_idx_by_rule[sort_by_idx]
        kmer_sequence_by_rule = kmer_sequence_by_rule[sort_by_idx]
        readdressed_kmer_idx_by_rule = dict((s, i) for i, s in enumerate(kmer_sequence_by_rule))
        readdressed_model = _readdress_model(model=model, rule_new_idx_by_kmer_seq=readdressed_kmer_idx_by_rule)
        X = np.vstack((_unpack_binary_bytes_from_ints(kmer_matrix[:, idx]) for idx in kmer_idx_by_rule)).T
        train_predictions = readdressed_model.predict(X[train_example_idx])
        progress_callback("Testing", 1.0 * len(train_example_idx) / (len(train_example_idx) + len(test_example_idx)))
        test_predictions = readdressed_model.predict(X[test_example_idx])
        progress_callback("Testing", 1.0)
    else:
        # The case where the model is just a leaf
        train_predictions = model.predict(np.empty((len(train_example_idx), 1)))
        progress_callback("Testing", 1.0 * len(train_example_idx) / (len(train_example_idx) + len(test_example_idx)))
        test_predictions = model.predict(np.empty((len(test_example_idx), 1)))
        progress_callback("Testing", 1.0)
    return train_predictions, test_predictions
    
def _learn_pruned_tree(hps, dataset_file, split_name):
    """
    Learns a cost-complexity pruned decision tree for a fixed set of hyperparameters and returns an estimate of its
    generalization error

    Parameters:
    -----------
    dataset_file: str
        The path to the Kover dataset
    split_name: str
        The name of the train/test split to use
    hps: dict
        A dictionnary of hyperparameter values (one value per key)

    Returns:
    --------
    cv_score: float
        The cross-validation score of the best tree that was grown using the provided HPs
    best_tree: float
        A tree that was grown from the entire training set. It is the one with the best generalization performance
        estimated by CV. It has also been pruned to prevent overfitting.

    Notes:
    ------
    The pruning is done using the cross-validation version of the cost-complexity pruning procedure described in:
    Breiman, L., Friedman, J., Stone, C. J., & Olshen, R. A. (1984). Classification and regression trees. CRC press.

    """

    # Open the dataset and load some split info into memory
    logging.debug("Opening the Kover dataset and loading split information into memory")
    dataset = KoverDataset(dataset_file)
    split = dataset.get_split(split_name)
    split.train_genome_idx = split.train_genome_idx[...]
    split.test_genome_idx = split.test_genome_idx[...]
    nb_classes = dataset.phenotype_tags.shape[0]


    # Load some dataset information into memory
    logging.debug("Loading dataset information into memory")
    rules = LazyKmerRuleList(dataset.kmer_sequences, dataset.kmer_by_matrix_column)
    rule_classifications = KmerRuleClassifications(dataset.kmer_matrix, dataset.genome_count)
    example_labels = dataset.phenotype.metadata[...]

    # Initialize the trees to be grown
    logging.debug("Planting seeds")
    fold_predictors = [DecisionTreeClassifier(criterion=hps["criterion"], max_depth=hps["max_depth"],
                                              min_samples_split=hps["min_samples_split"],
                                              class_importance=hps["class_importance"])
                       for _ in xrange(len(split.folds))]
    master_predictor = DecisionTreeClassifier(criterion=hps["criterion"], max_depth=hps["max_depth"],
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
                               example_idx = {c:fold.train_genome_idx[example_labels[fold.train_genome_idx] == c] for c in range(nb_classes)},
                               rule_blacklist=None,  # TODO: blacklist not implemented
                               tiebreaker=None,
                               level_callback=None,
                               split_callback=None)

    # Also build an overgrown decision tree on the entire dataset
    logging.debug("Growing the master tree")
    master_predictor.fit(rules=rules,
                         rule_classifications=rule_classifications,
                         example_idx = {c:split.train_genome_idx[example_labels[split.train_genome_idx] == c] for c in range(nb_classes)},
                         rule_blacklist=None,  # TODO: blacklist not implemented
                         tiebreaker=None,
                         level_callback=None,
                         split_callback=None)

    # Get the pruned master and cross-validation trees
    master_alphas, master_pruned_trees = _prune_tree(master_predictor.model.decision_tree)
    fold_alphas = []
    fold_pruned_trees = []
    for i in xrange(len(split.folds)):
		alphas, trees = _prune_tree(fold_predictors[i].model.decision_tree)
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
            fold_test_risk = _get_metrics(_predictions(t, dataset.kmer_matrix, [], fold_test_example_idx)[1],
                                                fold_example_labels)["risk"][0]
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

    for i in xrange(len(master_alphas) - 1):
        alpha_k = master_alphas[i]
        alpha_kplus1 = master_alphas[i + 1]
        alphaprime_k = sqrt(alpha_k * alpha_kplus1)
        cv_score = np.mean([fold_scores_by_alpha[j][alphaprime_k] for j in xrange(len(split.folds))])
        if cv_score <= min_score:  # Note: assumes that alphas are sorted in increasing order
            min_score = cv_score
            min_score_tree = master_pruned_trees[i]

    # Return the best tree and its error estimate
    return hps, min_score, min_score_tree
    
def learn(dataset_file, split_name, criterion, max_depth, min_samples_split, class_importance, no_tiebreaker,
          parameter_selection, n_cpu, progress_callback=None, warning_callback=None, error_callback=None):
    """


    """
    # Initialize callback functions
    warning_callback, error_callback, progress_callback = _init_callback_functions(warning_callback, error_callback,
                                                                                   progress_callback)
    # Check and initialize (hyper)parameters
    if n_cpu is None:
        n_cpu = cpu_count()
    criterion = np.unique(criterion)
    class_importance = np.unique(class_importance)
    max_depth = np.unique(max_depth)
    min_samples_split = np.unique(min_samples_split)

    # Find the best combination of hyperparameters
    pool = Pool(n_cpu)
    _hp_eval_func = partial(_learn_pruned_tree, dataset_file=dataset_file, split_name=split_name)
    best_hps = None
    best_score = np.infty
    best_tree = None
    for hps, score, tree in pool.imap_unordered(_hp_eval_func,
                                                ({"criterion": hps[0],
                                                  "class_importance": hps[1],
                                                  "max_depth": hps[2],
                                                  "min_samples_split": hps[3]}
                                                  for hps in product(criterion, class_importance, max_depth,
                                                                     min_samples_split))):
        # TODO: Add more logic, like prefering balanced class importances, etc.
        if score < best_score:
            best_hps = hps
            best_score = score
            best_tree = tree

    # Open the dataset and load some split info into memory
    logging.debug("Opening the Kover dataset and loading split information into memory")
    dataset = KoverDataset(dataset_file)
    split = dataset.get_split(split_name)
    split.train_genome_idx = split.train_genome_idx[...]
    split.test_genome_idx = split.test_genome_idx[...]
    example_labels = dataset.phenotype.metadata[...]

    # Using the best hyperparameters, compute predictions and metrics
    train_predictions, test_predictions = _predictions(best_tree, dataset.kmer_matrix, split.train_genome_idx,
                                                       split.test_genome_idx, progress_callback)

    train_answers = example_labels[split.train_genome_idx]
    test_answers = example_labels[split.test_genome_idx]

    train_metrics = _get_metrics(train_predictions, train_answers)
    if len(split.test_genome_idx) > 0:
        test_metrics = _get_metrics(test_predictions, test_answers)
    else:
        test_metrics = None

    # Get the idx of the training/testing examples that are correctly/incorrectly classified by the model
    classifications = defaultdict(list)
    classifications["train_correct"] = dataset.genome_identifiers[split.train_genome_idx[train_predictions == train_answers].tolist()].tolist() if train_metrics["risk"][0] < 1.0 else []
    classifications["train_errors"] = dataset.genome_identifiers[split.train_genome_idx[train_predictions != train_answers].tolist()].tolist() if train_metrics["risk"][0] > 0 else []
    if len(split.test_genome_idx) > 0:
        classifications["test_correct"] = dataset.genome_identifiers[split.test_genome_idx[test_predictions == test_answers].tolist()].tolist() if test_metrics["risk"][0] < 1.0 else []
        classifications["test_errors"] = dataset.genome_identifiers[split.test_genome_idx[test_predictions != test_answers].tolist()].tolist() if test_metrics["risk"][0] > 0 else []

    return best_hps, best_score, train_metrics, test_metrics, best_tree, dict((str(r), 1.0) for r in best_tree.rules), classifications
    # TODO: missing model_equivalent_rules, rule importances
