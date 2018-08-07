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
import numpy as np

from collections import defaultdict

def _get_binary_metrics(predictions, answers):
    """
    Compute metrics for binary predictions.
    Parameters:
    -----------
    predictions: numpy_array, shape=(n_predictions, n_examples)
	    The predictions array.
    answers: numpy_array, shape=(n_examples,)
	    The labels vector associated to the examples.
    Returns:
    --------
    metrics: dictionary
	    A dictionary containing for each metric key,
	    a vector of values associated to the predictions.
    """
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

def _get_multiclass_metrics(predictions, answers, nb_class):
    """
    Compute metrics for multi-class predictions.
    Parameters:
    -----------
    predictions: numpy_array, shape=(n_predictions, n_examples)
	    The predictions array.
    answers: numpy_array, shape=(n_examples,)
	    The labels vector associated to the examples.
    Returns:
    --------
    metrics: dictionary
	    A dictionary containing for each metric key,
	    a vector of values associated to the predictions.
    """
    if len(predictions.shape) == 1:
	    predictions = predictions.reshape(1, -1)
    metrics = defaultdict(list)
    for i in xrange(predictions.shape[0]):
	    p = predictions[i]
	    risk = 1.0 * len(p[p != answers]) / len(answers)
	    confusion_matrix = []
	    for actual_class in range(nb_class):
		    confusion_matrix.append([len(np.where(p[answers == actual_class] == predicted_class)[0]) \
										for predicted_class in range(nb_class)])
	    metrics["risk"].append(risk)
	    metrics["confusion_matrix"].append(confusion_matrix)
    return metrics
