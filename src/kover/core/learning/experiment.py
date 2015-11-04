__author__ = 'Alexandre'

import h5py as h
import logging
import numpy as np

from collections import defaultdict
from copy import deepcopy
from functools import partial
from itertools import product
from multiprocessing import Pool, cpu_count

from ..dataset.ds import KoverDataset
from .set_covering_machine.models import ConjunctionModel, DisjunctionModel
from .set_covering_machine.rules import LazyKmerRuleList, KmerRuleClassifications
from .set_covering_machine.scm import SetCoveringMachine
from ..utils import _duplicate_last_element, _unpack_binary_bytes_from_ints

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

def _predictions(model, kmer_matrix, train_example_idx, test_example_idx, progress_callback=None):
    """
    Makes predictions by loading only the columns of the kmer matrix that are targetted by the model.
    """
    if progress_callback is None:
        progress_callback = lambda t, p: None

    progress_callback("Testing", 0.0)
    kmer_idx_by_rule = np.array([r.kmer_index for r in model])
    kmer_sequence_by_rule = np.array([r.kmer_sequence for r in model])

    sort_by_idx = np.argsort(kmer_idx_by_rule)
    kmer_idx_by_rule = kmer_idx_by_rule[sort_by_idx]
    kmer_sequence_by_rule = kmer_sequence_by_rule[sort_by_idx]

    readdressed_kmer_idx_by_rule = dict((s, i) for i, s in enumerate(kmer_sequence_by_rule))
    readdressed_model = _readdress_model(model=model, kmer_idx=readdressed_kmer_idx_by_rule)

    X = _unpack_binary_bytes_from_ints(kmer_matrix[:, kmer_idx_by_rule.tolist()])
    train_predictions = readdressed_model.predict(X[train_example_idx])
    progress_callback("Testing", 1.0 * len(train_example_idx) / (len(train_example_idx) + len(test_example_idx)))
    test_predictions = readdressed_model.predict(X[test_example_idx])
    progress_callback("Testing", 1.0)

    return train_predictions, test_predictions

def _readdress_model(model, kmer_idx):
    """
    Produces a new model that looks for the k-mers at different locations in the input
    """
    new_model = deepcopy(model)
    for r in new_model:
        r.kmer_index = kmer_idx[r.kmer_sequence]
    return new_model

def print_hp(mt, p, max_rules):
    print "HP:", mt, p, max_rules
    from time import sleep
    np.random.seed()
    t = np.random.randint(0, 10)
    sleep(t)
    return (mt, p), 1, t

def _cv_score_hp(hp_values, max_rules, dataset_file, split_name):
    model_type = hp_values[0]
    p = hp_values[1]

    dataset = KoverDataset(dataset_file)
    folds = dataset.get_split(split_name).folds
    rules = LazyKmerRuleList(dataset.kmer_sequences, dataset.kmer_by_matrix_column)
    rule_classifications = KmerRuleClassifications(dataset.kmer_matrix, dataset.genome_count)

    def _iteration_callback(iteration_infos, tmp_model, test_predictions_by_model_length, test_example_idx):
        selected_rule = rules[iteration_infos["selected_attribute_idx"]]
        tmp_model.add(selected_rule)
        _, test_predictions = _predictions(tmp_model, dataset.kmer_matrix, [], test_example_idx)
        test_predictions_by_model_length.append(test_predictions)

    def _tiebreaker(best_utility_idx, attribute_classifications, positive_error_counts, negative_cover_counts,
               positive_example_idx, negative_example_idx, rule_risks, model_type):
        logging.debug("There are %d candidate attributes." % len(best_utility_idx))
        tie_rule_risks = rule_risks[best_utility_idx]
        if model_type == "conjunction":
            result = best_utility_idx[tie_rule_risks == tie_rule_risks.min()]
        else:
            # Use max instead of min, since in the disjunction case the risks = 1.0 - conjunction risks (inverted ys)
            result = best_utility_idx[tie_rule_risks == tie_rule_risks.max()]
        return result

    fold_score_by_model_length = np.zeros((len(folds), max_rules))
    for i, fold in enumerate(folds):
        logging.debug("Fold: %s" % fold.name)
        rule_risks = np.hstack((fold.unique_risk_by_kmer[...],
                                fold.unique_risk_by_anti_kmer[...])) # Too bad that we need to load each time. Maybe invert the loops (all hp for each fold)

        train_example_idx = fold.train_genome_idx
        test_example_idx = fold.test_genome_idx
        positive_example_idx = train_example_idx[dataset.phenotype.metadata[train_example_idx] == 1].reshape(-1)
        negative_example_idx = train_example_idx[dataset.phenotype.metadata[train_example_idx] == 0].reshape(-1)
        tiebreaker = partial(_tiebreaker, rule_risks=rule_risks, model_type=model_type)
        test_predictions_by_model_length = []
        tmp_model = ConjunctionModel() if model_type == "conjunction" else DisjunctionModel()
        iteration_callback = partial(_iteration_callback, tmp_model=tmp_model,
                                     test_predictions_by_model_length=test_predictions_by_model_length,
                                     test_example_idx=test_example_idx)

        predictor = SetCoveringMachine(model_type=model_type, p=p, max_rules=max_rules)
        predictor.fit(rules=rules,
                      rule_classifications=rule_classifications,
                      positive_example_idx=positive_example_idx,
                      negative_example_idx=negative_example_idx,
                      tiebreaker=tiebreaker,
                      iteration_callback=iteration_callback)
        test_predictions_by_model_length = np.array(_duplicate_last_element(test_predictions_by_model_length, max_rules))
        fold_score_by_model_length[i] = _get_metrics(test_predictions_by_model_length,
                                                       dataset.phenotype.metadata[test_example_idx])["risk"]

    score_by_model_length = np.mean(fold_score_by_model_length, axis=0)
    best_score_idx = np.argmin(score_by_model_length)
    best_hp_score = score_by_model_length[best_score_idx]
    best_model_length = best_score_idx + 1
    logging.debug("Returning")
    return (model_type, p, best_model_length), best_hp_score

def _cross_validation(dataset_file, split_name, model_types, p_values, max_rules, n_cpu, progress_callback,
                      warning_callback, error_callback):
    """
    Returns the best parameter combination and its cv score
    """
    n_hp_combinations = len(model_types) * len(p_values)
    logging.debug("There are %d hyperparameter combinations to try." % n_hp_combinations)

    logging.debug("Using %d CPUs." % n_cpu)
    pool = Pool(processes=n_cpu)
    hp_eval_func = partial(_cv_score_hp, dataset_file=dataset_file, split_name=split_name, max_rules=max_rules)

    best_hp_score = 1.0
    best_hp = {"model_type": None, "p": None, "max_rules": None}
    n_completed = 0.0
    progress_callback("Cross-validation", 0.0)
    for hp, score in pool.imap_unordered(hp_eval_func, product(model_types, p_values)):
        n_completed += 1
        progress_callback("Cross-validation", n_completed / n_hp_combinations)
        if (score < best_hp_score) or \
           (score == best_hp_score and hp[2] < best_hp["max_rules"]) or \
           (score == best_hp_score and abs(1.0 - hp[1]) < abs(1.0 - best_hp["p"])):
            best_hp["model_type"] = hp[0]
            best_hp["p"] = hp[1]
            best_hp["max_rules"] = hp[2]
            best_hp_score = score
    return best_hp_score, best_hp

def learn(dataset_file, split_name, model_type, p, max_rules, parameter_selection, n_cpu, progress_callback=None,
          warning_callback=None, error_callback=None):
    """
    parameter_selection: bound, cv, none (use first value of each if multiple)
    """
    # Execution callback functions
    if warning_callback is None:
        warning_callback = lambda w: logging.warning(w)
    if error_callback is None:
        def normal_raise(exception):
            raise exception
        error_callback = normal_raise
    if progress_callback is None:
        progress_callback = lambda t, p: None

    if n_cpu is None:
        n_cpu = cpu_count()

    model_type = np.unique(model_type)
    p = np.unique(p)

    dataset = KoverDataset(dataset_file)

    # Score the hyperparameter combinations
    # ------------------------------------------------------------------------------------------------------------------
    if parameter_selection == "bound":
        error_callback(NotImplementedError("Not supported yet!")) #TODO: Determine if we need the max genome size and delta.
    elif parameter_selection == "cv":
        n_folds = len(dataset.get_split(split_name).folds)
        if n_folds < 1:
            error_callback(Exception("Cross-validation cannot be performed on a split with no folds."))
        best_hp_score, best_hp = _cross_validation(dataset_file, split_name, model_type, p, max_rules, n_cpu,
                                                   progress_callback, warning_callback, error_callback)
    else:
        # Use the first value provided for each parameter
        best_hp = {"model_type": model_type[0], "p": p[0], "max_rules": max_rules}
        best_hp_score = None

    # Use the best hyperparameters to train/test on the split
    # ------------------------------------------------------------------------------------------------------------------
    full_train_progress = {"n_rules": 0.0}
    def _iteration_callback(iteration_infos):
        full_train_progress["n_rules"] += 1

    def _tiebreaker(best_utility_idx, attribute_classifications, positive_error_counts, negative_cover_counts,
           positive_example_idx, negative_example_idx, rule_risks, model_type):
        logging.debug("There are %d candidate attributes." % len(best_utility_idx))
        tie_rule_risks = rule_risks[best_utility_idx]
        if model_type == "conjunction":
            result = best_utility_idx[tie_rule_risks == tie_rule_risks.min()]
        else:
            # Use max instead of min, since in the disjunction case the risks = 1.0 - conjunction risks (inverted ys)
            result = best_utility_idx[tie_rule_risks == tie_rule_risks.max()]
        return result

    split = dataset.get_split(split_name)
    rules = LazyKmerRuleList(dataset.kmer_sequences, dataset.kmer_by_matrix_column)
    rule_classifications = KmerRuleClassifications(dataset.kmer_matrix, dataset.genome_count)

    train_example_idx = split.train_genome_idx
    test_example_idx = split.test_genome_idx
    positive_example_idx = train_example_idx[dataset.phenotype.metadata[train_example_idx] == 1].reshape(-1)
    negative_example_idx = train_example_idx[dataset.phenotype.metadata[train_example_idx] == 0].reshape(-1)

    predictor = SetCoveringMachine(model_type=best_hp["model_type"], p=best_hp["p"], max_rules=best_hp["max_rules"])
    progress_callback("Training", full_train_progress["n_rules"] / best_hp["max_rules"])
    predictor.fit(rules=rules,
                  rule_classifications=rule_classifications,
                  positive_example_idx=positive_example_idx,
                  negative_example_idx=negative_example_idx,
                  tiebreaker= partial(_tiebreaker,
                                      rule_risks=np.hstack((split.unique_risk_by_kmer[...],
                                                            split.unique_risk_by_anti_kmer[...])),
                                      model_type=model_type))

    train_predictions, test_predictions = _predictions(predictor.model, dataset.kmer_matrix, train_example_idx,
                                                       test_example_idx, progress_callback)
    train_metrics = _get_metrics(train_predictions, dataset.phenotype.metadata[train_example_idx])
    test_metrics = _get_metrics(test_predictions, dataset.phenotype.metadata[test_example_idx])

    return best_hp, best_hp_score, train_metrics, test_metrics, predictor.model