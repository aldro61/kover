#!/usr/bin/env python
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

import logging
import numpy as np

from collections import defaultdict
from copy import deepcopy
from functools import partial
from itertools import product
from math import exp, log as ln, pi
from multiprocessing import Pool, cpu_count
from scipy.misc import comb

from ...dataset.ds import KoverDataset
from ..common.models import ConjunctionModel, DisjunctionModel
from ..common.rules import LazyKmerRuleList, KmerRuleClassifications
from ..learners.scm import SetCoveringMachine
from ...utils import (
    _duplicate_last_element,
    _unpack_binary_bytes_from_ints,
    _parse_kmer_blacklist,
)
from ..experiments.metrics import _get_binary_metrics


def _predictions(
    model, kmer_matrix, train_example_idx, test_example_idx, progress_callback=None
):
    """Computes predictions by loading only the columns of the kmer matrix that are targetted by the model.

    Parameters
    ----------
    model: BaseModel
        The model used for predicting.

    kmer_matrix: BaseRuleClassifications
        The matrix containing the classifications of each rule on each learning example.

    train_example_idx: array-like, dtype=uint
        The index of the rows of kmer_matrix corresponding to the training examples.

    test_example_idx: array-like, dtype=uint
        The index of the rows of kmer_matrix corresponding to the testing examples.

    progress_callback: function with arguments task, percent_completed
        A callback function used to keep track of the task's completion.

    """
    if progress_callback is None:
        progress_callback = lambda t, p: None

    progress_callback("Testing", 0.0)

    if len(model) == 0:
        train_predictions = model.predict(np.zeros(len(train_example_idx)))
        test_predictions = model.predict(np.zeros(len(test_example_idx)))
    else:
        # We use h5py to load only the columns of the k-mer matrix targeted by the model. The indices passed to h5py
        # need to be sorted. We change the kmer_idx of the rules in the model to be 0 ... n_rules, with the rule that
        # initially had the smallest kmer_idx pointing to 0 and the one with the largest kmer_idx pointing to n_rules.
        # We then load only the appropriate columns and apply the readdressed model to the data (in RAM).
        columns_to_load = []
        readdressed_model = deepcopy(model)
        for i, rule_idx in enumerate(np.argsort([r.kmer_index for r in model.rules])):
            rule = readdressed_model.rules[rule_idx]
            columns_to_load.append(rule.kmer_index)
            rule.kmer_index = i

        # Load the columns targeted by the model and make predictions using the readdressed model
        X = _unpack_binary_bytes_from_ints(kmer_matrix[:, columns_to_load])
        train_predictions = readdressed_model.predict(X[train_example_idx])
        progress_callback(
            "Testing",
            1.0
            * len(train_example_idx)
            / (len(train_example_idx) + len(test_example_idx)),
        )
        test_predictions = readdressed_model.predict(X[test_example_idx])

    progress_callback("Testing", 1.0)

    return train_predictions, test_predictions


def _cv_score_hp(hp_values, max_rules, dataset_file, split_name, rule_blacklist):
    model_type = hp_values[0]
    p = hp_values[1]

    dataset = KoverDataset(dataset_file)
    folds = dataset.get_split(split_name).folds
    rules = LazyKmerRuleList(dataset.kmer_sequences, dataset.kmer_by_matrix_column)
    rule_classifications = KmerRuleClassifications(
        dataset.kmer_matrix, dataset.genome_count
    )

    def _iteration_callback(
        iteration_infos, tmp_model, test_predictions_by_model_length, test_example_idx
    ):
        tmp_model.add(iteration_infos["selected_rule"])
        _, test_predictions = _predictions(
            tmp_model, dataset.kmer_matrix, [], test_example_idx
        )
        test_predictions_by_model_length.append(test_predictions)

    def _tiebreaker(best_utility_idx, rule_risks, model_type):
        logging.debug("There are %d candidate rules." % len(best_utility_idx))
        tie_rule_risks = rule_risks[best_utility_idx]
        if model_type == "conjunction":
            result = best_utility_idx[np.isclose(tie_rule_risks, tie_rule_risks.min())]
        else:
            # Use max instead of min, since in the disjunction case the risks = 1.0 - conjunction risks (inverted ys)
            result = best_utility_idx[np.isclose(tie_rule_risks, tie_rule_risks.max())]
        return result

    fold_score_by_model_length = np.ones((len(folds), max_rules + 1)) * np.infty
    for i, fold in enumerate(folds):
        logging.debug("Fold: %s" % fold.name)
        rule_risks = np.hstack(
            (fold.unique_risk_by_kmer[...], fold.unique_risk_by_anti_kmer[...])
        )  # Too bad that we need to load each time. Maybe invert the loops (all hp for each fold)

        train_example_idx = fold.train_genome_idx
        test_example_idx = fold.test_genome_idx
        positive_example_idx = train_example_idx[
            dataset.phenotype.metadata[train_example_idx] == 1
        ].reshape(-1)
        negative_example_idx = train_example_idx[
            dataset.phenotype.metadata[train_example_idx] == 0
        ].reshape(-1)
        tiebreaker = partial(_tiebreaker, rule_risks=rule_risks, model_type=model_type)
        test_predictions_by_model_length = []
        tmp_model = (
            ConjunctionModel() if model_type == "conjunction" else DisjunctionModel()
        )
        iteration_callback = partial(
            _iteration_callback,
            tmp_model=tmp_model,
            test_predictions_by_model_length=test_predictions_by_model_length,
            test_example_idx=test_example_idx,
        )

        predictor = SetCoveringMachine(model_type=model_type, p=p, max_rules=max_rules)

        # Empty model predictions (length = 0)
        # Do this before fitting in case there are no predictive features in the data
        test_predictions_by_model_length.append(
            _predictions(tmp_model, dataset.kmer_matrix, [], test_example_idx)[1]
        )

        predictor.fit(
            rules=rules,
            rule_classifications=rule_classifications,
            positive_example_idx=positive_example_idx,
            negative_example_idx=negative_example_idx,
            rule_blacklist=rule_blacklist,
            tiebreaker=tiebreaker,
            iteration_callback=iteration_callback,
        )

        # Calcule the risk for each model length
        # Note: If the model stopped adding rules before the max, then we use the score
        #       for the last added rule as the score for all subsequent lengths.
        test_predictions_by_model_length = np.array(
            _duplicate_last_element(test_predictions_by_model_length, max_rules + 1)
        )
        fold_score_by_model_length[i] = _get_binary_metrics(
            predictions=test_predictions_by_model_length,
            answers=dataset.phenotype.metadata[test_example_idx],
        )["risk"]

    score_by_model_length = np.mean(fold_score_by_model_length, axis=0)
    best_score_idx = np.argmin(score_by_model_length)
    best_hp_score = score_by_model_length[best_score_idx]
    best_model_length = best_score_idx

    return (model_type, p, best_model_length), best_hp_score


def _cross_validation(
    dataset_file,
    split_name,
    model_types,
    p_values,
    max_rules,
    rule_blacklist,
    n_cpu,
    progress_callback,
    warning_callback,
    error_callback,
):
    """
    Returns the best parameter combination and its cv score
    """
    n_hp_combinations = len(model_types) * len(p_values)
    logging.debug(
        "There are %d hyperparameter combinations to try." % n_hp_combinations
    )

    logging.debug("Using %d CPUs." % n_cpu)
    pool = Pool(processes=n_cpu)
    hp_eval_func = partial(
        _cv_score_hp,
        dataset_file=dataset_file,
        split_name=split_name,
        max_rules=max_rules,
        rule_blacklist=rule_blacklist,
    )

    best_hp_score = 1.0
    best_hp = {"model_type": None, "p": None, "max_rules": None}
    n_completed = 0.0
    progress_callback("Cross-validation", 0.0)
    for hp, score in pool.imap_unordered(hp_eval_func, product(model_types, p_values)):
        n_completed += 1
        progress_callback("Cross-validation", n_completed / n_hp_combinations)
        if (
            (not np.allclose(score, best_hp_score) and score < best_hp_score)
            or (np.allclose(score, best_hp_score) and hp[2] < best_hp["max_rules"])
            or (
                np.allclose(score, best_hp_score)
                and hp[2] == best_hp["max_rules"]
                and not np.allclose(hp[1], best_hp["p"])
                and abs(1.0 - hp[1]) < abs(1.0 - best_hp["p"])
            )
        ):
            best_hp["model_type"] = hp[0]
            best_hp["p"] = hp[1]
            best_hp["max_rules"] = hp[2]
            best_hp_score = score

    return best_hp_score, best_hp


def _full_train(
    dataset,
    split_name,
    model_type,
    p,
    max_rules,
    max_equiv_rules,
    rule_blacklist,
    random_generator,
    progress_callback,
):
    full_train_progress = {"n_rules": 0.0}

    def _iteration_callback(iteration_infos, model_type, equivalent_rules):
        full_train_progress["n_rules"] += 1
        progress_callback("Training", full_train_progress["n_rules"] / max_rules)

        # Ensure that there are no more equivalent rules than the specified maximum
        if len(iteration_infos["equivalent_rules_idx"]) > max_equiv_rules:
            logging.debug(
                "There are more equivalent rules than the allowed maximum. Subsampling %d rules."
                % max_equiv_rules
            )
            random_idx = random_generator.choice(
                len(iteration_infos["equivalent_rules_idx"]),
                max_equiv_rules,
                replace=False,
            )
            random_idx.sort()
            iteration_infos["equivalent_rules_idx"] = iteration_infos[
                "equivalent_rules_idx"
            ][random_idx]

        # Adjust and store the equivalent rule indices
        if model_type == "disjunction":
            n_kmers = rule_classifications.shape[1] / 2
            iteration_infos["equivalent_rules_idx"] += n_kmers
            iteration_infos["equivalent_rules_idx"] %= 2 * n_kmers
            equivalent_rules.append(iteration_infos["equivalent_rules_idx"])
        else:
            equivalent_rules.append(iteration_infos["equivalent_rules_idx"])

    def _tiebreaker(best_utility_idx, rule_risks, model_type):
        logging.debug("There are %d candidate rules." % len(best_utility_idx))
        tie_rule_risks = rule_risks[best_utility_idx]
        if model_type == "conjunction":
            result = best_utility_idx[np.isclose(tie_rule_risks, tie_rule_risks.min())]
        else:
            # Use max instead of min, since in the disjunction case the risks = 1.0 - conjunction risks (inverted ys)
            result = best_utility_idx[np.isclose(tie_rule_risks, tie_rule_risks.max())]
        return result

    rules = LazyKmerRuleList(dataset.kmer_sequences, dataset.kmer_by_matrix_column)
    rule_classifications = KmerRuleClassifications(
        dataset.kmer_matrix, dataset.genome_count
    )
    split = dataset.get_split(split_name)

    train_example_idx = split.train_genome_idx
    positive_example_idx = train_example_idx[
        dataset.phenotype.metadata[train_example_idx] == 1
    ].reshape(-1)
    negative_example_idx = train_example_idx[
        dataset.phenotype.metadata[train_example_idx] == 0
    ].reshape(-1)

    model_equivalent_rules = []
    predictor = SetCoveringMachine(model_type=model_type, p=p, max_rules=max_rules)

    # In some cases, the CV will find that an empty model works better.
    # There is no need to train in that case. We return dummy values.
    if max_rules == 0:
        return predictor.model, np.array([]), np.array([])

    progress_callback("Training", 0)
    predictor.fit(
        rules=rules,
        rule_classifications=rule_classifications,
        positive_example_idx=positive_example_idx,
        negative_example_idx=negative_example_idx,
        rule_blacklist=rule_blacklist,
        tiebreaker=partial(
            _tiebreaker,
            rule_risks=np.hstack(
                (split.unique_risk_by_kmer[...], split.unique_risk_by_anti_kmer[...])
            ),
            model_type=model_type,
        ),
        iteration_callback=partial(
            _iteration_callback,
            model_type=model_type,
            equivalent_rules=model_equivalent_rules,
        ),
    )

    return predictor.model, predictor.rule_importances, model_equivalent_rules


def _bound(
    train_predictions,
    train_answers,
    train_example_idx,
    model,
    delta,
    max_genome_size,
    rule_classifications,
):
    # Construct the smallest possible compression set (Chvatal greedy approx for minimum set cover)
    logging.debug("Constructing the compression set.")
    compression_set = []
    if len(model) > 0:
        presence_by_example = rule_classifications.get_columns(
            [r.kmer_index for r in model]
        )[train_example_idx]
        while presence_by_example.shape[1] != 0:
            score = presence_by_example.sum(axis=1)
            best_example_relative_idx = np.argmax(score)
            compression_set.append(best_example_relative_idx)
            presence_by_example = presence_by_example[
                :, presence_by_example[best_example_relative_idx] == 0
            ]
    logging.debug("The compression set contains %d examples." % len(compression_set))

    # Compute the bound value
    logging.debug("Computing the bound value.")
    h_card = float(len(model))
    Z_card = float(len(compression_set) * max_genome_size)
    m = float(len(train_answers))
    mz = float(len(compression_set))
    r = float(
        (train_predictions != train_answers).sum()
        - (train_predictions[compression_set] != train_answers[compression_set]).sum()
    )
    return 1.0 - exp(
        (-1.0 / (m - mz - r))
        * (
            ln(comb(m, mz, exact=True)) + ln(comb(m - mz, r, exact=True)) + 0
            if h_card == 0
            else (h_card * ln(2 * Z_card))
            + ln(
                pi ** 6
                * (h_card + 1) ** 2
                * (r + 1) ** 2
                * (mz + 1) ** 2
                / (216 * delta)
            )
        )
    )


def _bound_score_hp(
    hp_values,
    max_rules,
    dataset_file,
    split_name,
    max_equiv_rules,
    rule_blacklist,
    bound_delta,
    bound_max_genome_size,
    random_generator,
):
    model_type = hp_values[0]
    p = hp_values[1]

    dataset = KoverDataset(dataset_file)
    rules = LazyKmerRuleList(dataset.kmer_sequences, dataset.kmer_by_matrix_column)
    rule_classifications = KmerRuleClassifications(
        dataset.kmer_matrix, dataset.genome_count
    )

    def _iteration_callback(
        iteration_infos,
        tmp_model,
        train_example_idx,
        train_answers,
        score_by_length,
        model_by_length,
        equivalent_rules,
        rule_importances,
        rule_classifications,
    ):
        tmp_model.add(iteration_infos["selected_rule"])
        model_by_length.append(deepcopy(tmp_model))
        rule_importances.append(iteration_infos["rule_importances"])

        # Store equivalent rules
        # Ensure that there are no more equivalent rules than the specified maximum
        if len(iteration_infos["equivalent_rules_idx"]) > max_equiv_rules:
            logging.debug(
                "There are more equivalent rules than the allowed maximum. Subsampling %d rules."
                % max_equiv_rules
            )
            random_idx = random_generator.choice(
                len(iteration_infos["equivalent_rules_idx"]),
                max_equiv_rules,
                replace=False,
            )
            random_idx.sort()
            iteration_infos["equivalent_rules_idx"] = iteration_infos[
                "equivalent_rules_idx"
            ][random_idx]

        # Adjust and store the equivalent rule indices
        if model_type == "disjunction":
            n_kmers = rule_classifications.shape[1] / 2
            iteration_infos["equivalent_rules_idx"] += n_kmers
            iteration_infos["equivalent_rules_idx"] %= 2 * n_kmers
            equivalent_rules.append(iteration_infos["equivalent_rules_idx"])
        else:
            equivalent_rules.append(iteration_infos["equivalent_rules_idx"])

        # Compute the bound value for the current model length
        _, train_predictions = _predictions(
            tmp_model, dataset.kmer_matrix, [], train_example_idx
        )
        score_by_length[iteration_infos["iteration_number"] - 1] = _bound(
            train_predictions=train_predictions,
            train_answers=train_answers,
            train_example_idx=train_example_idx,
            model=tmp_model,
            delta=bound_delta,
            max_genome_size=bound_max_genome_size,
            rule_classifications=rule_classifications,
        )

    def _tiebreaker(best_utility_idx, rule_risks, model_type):
        logging.debug("There are %d candidate rules." % len(best_utility_idx))
        tie_rule_risks = rule_risks[best_utility_idx]
        if model_type == "conjunction":
            result = best_utility_idx[np.isclose(tie_rule_risks, tie_rule_risks.min())]
        else:
            # Use max instead of min, since in the disjunction case the risks = 1.0 - conjunction risks (inverted ys)
            result = best_utility_idx[np.isclose(tie_rule_risks, tie_rule_risks.max())]
        return result

    split = dataset.get_split(split_name)
    rule_risks = np.hstack(
        (split.unique_risk_by_kmer[...], split.unique_risk_by_anti_kmer[...])
    )
    train_example_idx = split.train_genome_idx
    positive_example_idx = train_example_idx[
        dataset.phenotype.metadata[train_example_idx] == 1
    ].reshape(-1)
    negative_example_idx = train_example_idx[
        dataset.phenotype.metadata[train_example_idx] == 0
    ].reshape(-1)
    train_answers = dataset.phenotype.metadata[train_example_idx]

    tiebreaker = partial(_tiebreaker, rule_risks=rule_risks, model_type=model_type)

    tmp_model = (
        ConjunctionModel() if model_type == "conjunction" else DisjunctionModel()
    )
    score_by_length = np.ones(max_rules)
    model_by_length = []
    equivalent_rules = []
    rule_importances = []
    iteration_callback = partial(
        _iteration_callback,
        tmp_model=tmp_model,
        train_example_idx=train_example_idx,
        train_answers=train_answers,
        score_by_length=score_by_length,
        model_by_length=model_by_length,
        equivalent_rules=equivalent_rules,
        rule_importances=rule_importances,
        rule_classifications=rule_classifications,
    )

    predictor = SetCoveringMachine(model_type=model_type, p=p, max_rules=max_rules)
    predictor.fit(
        rules=rules,
        rule_classifications=rule_classifications,
        positive_example_idx=positive_example_idx,
        negative_example_idx=negative_example_idx,
        rule_blacklist=rule_blacklist,
        tiebreaker=tiebreaker,
        iteration_callback=iteration_callback,
        iteration_rule_importances=True,
    )

    # Handle edge case where no rules were added to the model
    if len(tmp_model) == 0:
        _, train_predictions = _predictions(
            tmp_model, dataset.kmer_matrix, [], train_example_idx
        )
        bound_value = _bound(
            train_predictions=train_predictions,
            train_answers=train_answers,
            train_example_idx=train_example_idx,
            model=tmp_model,
            delta=bound_delta,
            max_genome_size=bound_max_genome_size,
            rule_classifications=rule_classifications,
        )
        best_model_length = 0
        best_hp_score = bound_value
        best_model = tmp_model
        best_rule_importances = np.array([])
        best_equivalent_rules = np.array([])
    else:
        best_score_idx = np.argmin(score_by_length)
        best_hp_score = score_by_length[best_score_idx]
        best_model = model_by_length[best_score_idx]
        best_rule_importances = rule_importances[best_score_idx]
        best_equivalent_rules = equivalent_rules[: best_score_idx + 1]
        best_model_length = best_score_idx + 1

    return (
        (model_type, p, best_model_length),
        best_hp_score,
        best_model,
        best_rule_importances,
        best_equivalent_rules,
    )


def _bound_selection(
    dataset_file,
    split_name,
    model_types,
    p_values,
    max_rules,
    max_equiv_rules,
    rule_blacklist,
    bound_delta,
    bound_max_genome_size,
    n_cpu,
    random_generator,
    progress_callback,
    warning_callback,
    error_callback,
):
    n_hp_combinations = len(model_types) * len(p_values)
    logging.debug(
        "There are %d hyperparameter combinations to try." % n_hp_combinations
    )

    logging.debug("Using %d CPUs." % n_cpu)
    pool = Pool(processes=n_cpu)
    hp_eval_func = partial(
        _bound_score_hp,
        dataset_file=dataset_file,
        split_name=split_name,
        max_rules=max_rules,
        max_equiv_rules=max_equiv_rules,
        rule_blacklist=rule_blacklist,
        bound_delta=bound_delta,
        bound_max_genome_size=bound_max_genome_size,
        random_generator=random_generator,
    )

    best_hp_score = 1.0
    best_hp = {"model_type": None, "p": None, "max_rules": None}
    n_completed = 0.0
    progress_callback("Bound selection", 0.0)
    for hp, score, model, rule_importances, equiv_rules in pool.imap_unordered(
        hp_eval_func, product(model_types, p_values)
    ):
        n_completed += 1
        progress_callback("Bound selection", n_completed / n_hp_combinations)
        if (
            (score < best_hp_score)
            or (score == best_hp_score and hp[2] < best_hp["max_rules"])
            or (
                score == best_hp_score
                and hp[2] == best_hp["max_rules"]
                and abs(1.0 - hp[1]) < abs(1.0 - best_hp["p"])
            )
        ):
            best_hp["model_type"] = hp[0]
            best_hp["p"] = hp[1]
            best_hp["max_rules"] = hp[2]
            best_hp_score = score
            best_model = model
            best_equiv_rules = equiv_rules
            best_rule_importances = rule_importances

    return best_hp_score, best_hp, best_model, best_rule_importances, best_equiv_rules


def _find_rule_blacklist(dataset_file, kmer_blacklist_file, warning_callback):
    """
    Finds the index of the rules that must be blacklisted.
    """
    dataset = KoverDataset(dataset_file)

    # Find all rules to blacklist
    rule_blacklist = []
    if kmer_blacklist_file is not None:
        kmers_to_blacklist = _parse_kmer_blacklist(
            kmer_blacklist_file, dataset.kmer_length
        )

        if kmers_to_blacklist:
            # XXX: the k-mers are assumed to be upper-cased in the dataset
            kmer_sequences = dataset.kmer_sequences[...].tolist()
            kmer_by_matrix_column = dataset.kmer_by_matrix_column[
                ...
            ].tolist()  # XXX: each k-mer is there only once (see wiki)
            n_kmers = len(kmer_sequences)

            kmers_not_found = []
            for k in kmers_to_blacklist:
                k = k.upper()
                try:
                    presence_rule_idx = kmer_by_matrix_column.index(
                        kmer_sequences.index(k)
                    )
                    absence_rule_idx = presence_rule_idx + n_kmers
                    rule_blacklist += [presence_rule_idx, absence_rule_idx]
                except ValueError:
                    kmers_not_found.append(k)

            if len(kmers_not_found) > 0:
                warning_callback(
                    "The following kmers could not be found in the dataset: "
                    + ", ".join(kmers_not_found)
                )

    return rule_blacklist


def learn_SCM(
    dataset_file,
    split_name,
    model_type,
    p,
    kmer_blacklist_file,
    max_rules,
    max_equiv_rules,
    parameter_selection,
    n_cpu,
    random_seed,
    authorized_rules,
    bound_delta=None,
    bound_max_genome_size=None,
    progress_callback=None,
    warning_callback=None,
    error_callback=None,
):
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

    random_generator = np.random.RandomState(random_seed)

    model_type = np.unique(model_type)
    p = np.unique(p)

    logging.debug("Searching for blacklisted rules.")
    rule_blacklist = _find_rule_blacklist(
        dataset_file=dataset_file,
        kmer_blacklist_file=kmer_blacklist_file,
        warning_callback=warning_callback,
    )

    dataset = KoverDataset(dataset_file)

    # Score the hyperparameter combinations
    # ------------------------------------------------------------------------------------------------------------------
    if parameter_selection == "bound":
        if bound_delta is None or bound_max_genome_size is None:
            error_callback(
                Exception(
                    "Bound selection cannot be performed without delta and the maximum genome length."
                )
            )

        # For bound selection, there is no need to retrain the algorithm after selecting the best hyperparameters.
        # The model is already obtained from all the training data. This is why we save the model here.
        (
            best_hp_score,
            best_hp,
            best_model,
            best_rule_importances,
            best_predictor_equiv_rules,
        ) = _bound_selection(
            dataset_file=dataset_file,
            split_name=split_name,
            model_types=model_type,
            p_values=p,
            max_rules=max_rules,
            rule_blacklist=rule_blacklist,
            max_equiv_rules=max_equiv_rules,
            bound_delta=bound_delta,
            bound_max_genome_size=bound_max_genome_size,
            n_cpu=n_cpu,
            random_generator=random_generator,
            progress_callback=progress_callback,
            warning_callback=warning_callback,
            error_callback=error_callback,
        )

    elif parameter_selection == "cv":
        n_folds = len(dataset.get_split(split_name).folds)
        if n_folds < 1:
            error_callback(
                Exception(
                    "Cross-validation cannot be performed on a split with no folds."
                )
            )
        best_hp_score, best_hp = _cross_validation(
            dataset_file=dataset_file,
            split_name=split_name,
            model_types=model_type,
            p_values=p,
            max_rules=max_rules,
            rule_blacklist=rule_blacklist,
            n_cpu=n_cpu,
            progress_callback=progress_callback,
            warning_callback=warning_callback,
            error_callback=error_callback,
        )

    else:
        # Use the first value provided for each parameter
        best_hp = {"model_type": model_type[0], "p": p[0], "max_rules": max_rules}
        best_hp_score = None

    # Use the best hyperparameters to train/test on the split
    # ------------------------------------------------------------------------------------------------------------------
    if parameter_selection == "bound":
        model = best_model
        equivalent_rules = best_predictor_equiv_rules
        rule_importances = best_rule_importances
    else:
        model, rule_importances, equivalent_rules = _full_train(
            dataset=dataset,
            split_name=split_name,
            model_type=best_hp["model_type"],
            p=best_hp["p"],
            max_rules=best_hp["max_rules"],
            max_equiv_rules=max_equiv_rules,
            rule_blacklist=rule_blacklist,
            random_generator=random_generator,
            progress_callback=progress_callback,
        )

    split = dataset.get_split(split_name)
    train_example_idx = split.train_genome_idx
    test_example_idx = split.test_genome_idx

    train_predictions, test_predictions = _predictions(
        model=model,
        kmer_matrix=dataset.kmer_matrix,
        train_example_idx=train_example_idx,
        test_example_idx=test_example_idx,
        progress_callback=progress_callback,
    )

    train_answers = dataset.phenotype.metadata[train_example_idx]
    train_metrics = _get_binary_metrics(train_predictions, train_answers)

    # No need to recompute the bound if bound selection was used
    if parameter_selection == "bound":
        train_metrics["bound"] = best_hp_score
    else:
        train_metrics["bound"] = _bound(
            train_predictions=train_predictions,
            train_answers=train_answers,
            train_example_idx=train_example_idx,
            model=model,
            delta=bound_delta,
            max_genome_size=bound_max_genome_size,
            rule_classifications=KmerRuleClassifications(
                dataset.kmer_matrix, dataset.genome_count
            ),
        )

    # Test metrics are computed only if there is a testing set
    if len(test_example_idx) > 0:
        test_answers = dataset.phenotype.metadata[test_example_idx]
        test_metrics = _get_binary_metrics(test_predictions, test_answers)
    else:
        test_metrics = None

    # Get the idx of the training/testing examples that are correctly/incorrectly classified by the model
    classifications = defaultdict(list)
    classifications["train_correct"] = (
        dataset.genome_identifiers[
            train_example_idx[train_predictions == train_answers].tolist()
        ].tolist()
        if train_metrics["risk"][0] < 1.0
        else []
    )
    classifications["train_errors"] = (
        dataset.genome_identifiers[
            train_example_idx[train_predictions != train_answers].tolist()
        ].tolist()
        if train_metrics["risk"][0] > 0
        else []
    )
    if len(test_example_idx) > 0:
        classifications["test_correct"] = (
            dataset.genome_identifiers[
                test_example_idx[test_predictions == test_answers].tolist()
            ].tolist()
            if test_metrics["risk"][0] < 1.0
            else []
        )
        classifications["test_errors"] = (
            dataset.genome_identifiers[
                test_example_idx[test_predictions != test_answers].tolist()
            ].tolist()
            if test_metrics["risk"][0] > 0
            else []
        )

    # Convert the equivalent rule indexes to rule objects
    rules = LazyKmerRuleList(dataset.kmer_sequences, dataset.kmer_by_matrix_column)
    model_equivalent_rules = [
        [rules[i] for i in equiv_idx] for equiv_idx in equivalent_rules
    ]

    return (
        best_hp,
        best_hp_score,
        train_metrics,
        test_metrics,
        model,
        rule_importances,
        model_equivalent_rules,
        classifications,
    )
