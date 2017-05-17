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

import numpy as np

conjunction = "conjunction"
disjunction = "disjunction"
scm = "scm"
cart = "cart"

class BaseModel(object):
    def __init__(self):
        self.rules = None
        super(BaseModel, self).__init__()
        
    def predict(self, X):
        raise NotImplementedError()
        
    @property
    def learner(self):
        raise NotImplementedError()
        
    def _to_string(self, separator=" "):
        raise NotImplementedError()
        
    def __str__(self):
        return self._to_string()
        
class CART_Model(BaseModel):
    def __init__(self):
        super(CART_Model, self).__init__()
    
    @property
    def learner(self):
        return cart

class SCM_Model(BaseModel):
    def __init__(self):
        super(SCM_Model, self).__init__()
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

class ConjunctionModel(SCM_Model):
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


class DisjunctionModel(SCM_Model):
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
