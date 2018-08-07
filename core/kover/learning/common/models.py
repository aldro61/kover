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

conjunction = "conjunction"
disjunction = "disjunction"
scm = "scm"
cart = "cart"

class BaseModel(object):
    def __init__(self):
        super(BaseModel, self).__init__()

    def predict(self, X):
        raise NotImplementedError()

	def predict_proba(self, X):
	    raise NotImplementedError()

    @property
    def learner(self):
        raise NotImplementedError()

    def __str__(self):
        return self._to_string()


class CARTModel(BaseModel):
    def __init__(self, class_tags=None):
        super(CARTModel, self).__init__()
        self.decision_tree = None
        self.class_tags = class_tags

    def predict(self, X):
        if self.decision_tree is None:
            raise RuntimeError("A decision tree must be fitted prior to calling predict.")
        predictions = self.decision_tree.predict(X)
        return np.asarray(predictions, dtype=np.uint8)

    def predict_proba(self, X):
        if self.decision_tree is None:
            raise RuntimeError("A decision tree must be fitted prior to calling predict.")
        self.decision_tree.predict_proba(X)

    @property
    def learner(self):
        return cart

    def _to_string(self, node=None,depth=0):
        if node is None:
            if self.decision_tree is None:
                print("No tree has been added to the model")
            node = self.decision_tree

        if self.class_tags is None:
            return str(self.decision_tree)

        tree_str = ""

        # Case : node is a leaf
        if node.is_leaf:
			tree_str += "\n" + ("    " * depth) + str(self.class_tags[node.class_prediction])

        # Case : node has two children
        else:
            # Print right branch
            tree_str += self._to_string(node=node.right_child, depth=depth + 1)

            # Print own value
            tree_str += "\n" + ("    " * depth + "   ") + str("/")
            tree_str += "\n" + ("    " * depth) + str(node.rule)
            tree_str += "\n" + ("    " * depth + "   ") + str("\\")

            # Print left_child branch
            tree_str += self._to_string(node=node.left_child, depth=depth + 1)

        return tree_str

    def __len__(self):
        if self.decision_tree is None:
            return 0
        return len(self.decision_tree)

    @property
    def depth(self):
        if self.decision_tree is None:
            return 0
        return self.decision_tree.tree_depth


class SCMModel(BaseModel):
    def __init__(self):
        super(SCMModel, self).__init__()
        self.rules = []

    def add(self, rule):
        self.rules.append(rule)

    def predict(self, X):
        predictions = self.predict_proba(X)
        predictions[predictions > 0.5] = 1
        predictions[predictions <= 0.5] = 0
        return np.asarray(predictions, dtype=np.uint8)

    def predict_proba(self, X):
        raise NotImplementedError()

    def remove(self, index):
        del self.rules[index]

    @property
    def example_dependencies(self):
        return [d for ba in self.rules for d in ba.example_dependencies]

    @property
    def learner(self):
        return scm

    @property
    def type(self):
        raise NotImplementedError()

    def _to_string(self, separator=" "):
        return separator.join([str(a) for a in self.rules])

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    def __iter__(self):
        for ba in self.rules:
            yield ba

    def __len__(self):
        return len(self.rules)


class ConjunctionModel(SCMModel):
    def predict_proba(self, X):
        predictions = np.ones(X.shape[0], np.float32)
        for a in self.rules:
            predictions *= a.classify(X)
        return predictions

    @property
    def type(self):
        return conjunction

    def __str__(self):
        return self._to_string(separator=" and ")


class DisjunctionModel(SCMModel):
    def predict_proba(self, X):
        predictions = np.ones(X.shape[0], dtype=np.float32)
        for a in self.rules:
            predictions *= 1.0 - a.classify(X) # Proportion of the voters that predict False in a
        return 1.0 - predictions

    @property
    def type(self):
        return disjunction

    def __str__(self):
        return self._to_string(separator=" or ")
