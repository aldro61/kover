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

from collections import defaultdict, deque 
from math import ceil 

UTIL_BLOCK_SIZE = 1000000

class BreimanInfo(object):     
	def __init__(self, node_n_examples_by_class, class_priors, total_n_examples_by_class):
		# Eq. 2.2 Probability that an example is in class j and falls into node t         
		self.p_j_t = [pi_j * N_j_t / N_j for pi_j, N_j_t, N_j in zip(class_priors, node_n_examples_by_class, 
																		  total_n_examples_by_class)] 
		# Eq. 2.3 Probability that any example falls in node t         
		self.p_t = sum(self.p_j_t)         
		# Eq. 2.4 Probability that an example is in class j given that it falls in node t         
		self.p_j_given_t = [p_j_t / self.p_t for p_j_t in self.p_j_t]         
		# Def. 2.10 Probability of misclassification given that an example falls into node t         
		self.r_t = 1.0 - max(self.p_j_given_t)         
		# Contribution of the node to the tree's overall missclassification rate         
		self.R_t = self.r_t * self.p_t
		
class TreeNode(object):
	# TODO: implement memoizing for attributes that require longer computation     
	def __init__(self, depth, class_examples_idx, total_n_examples_by_class, class_priors, rule=None, parent=None, left_child=None,                  
				right_child=None, criterion_value=-1):         
		self.rule = rule         
		self.parent = parent         
		self.left_child = left_child         
		self.right_child = right_child         
		self.class_examples_idx = class_examples_idx         
		self.depth = depth         
		self.criterion_value = criterion_value         
		self.breiman_info = BreimanInfo([len(class_examples_idx[c]) for c in range(len(class_examples_idx))],                                         
		class_priors, total_n_examples_by_class)
		
		
	@property     
	def is_leaf(self):         
		return self.left_child is None and self.right_child is None     
		
	@property     
	def is_root(self):         
		return self.parent is None and self.left_child is not None and self.right_child is not None
		
	@property     
	def n_examples(self):         
		"""        
		Returns the number of examples in the node         
		"""   
		return  sum(len(x) for x in self.class_examples_idx.itervalues())       
		
	@property     
	def class_proportions(self):         
		"""         
		Returns the proportion of examples of each class in the node         
		"""  
		return {key:1.0 * len(value)/self.n_examples for key, value in self.class_examples_idx.items()}
	
	@property     
	def rules(self):        
		return _get_tree_rules(self)
	
	def __iter__(self):         
		for r in _get_tree_rules(self):             
			yield r
	
	def __len__(self):         
		"""         
		Returns the number of rules in the tree         
		"""         
		return len(_get_tree_rules(self))
	
		
	def _get_tree_leaves(root):     
		def _get_leaves(node):         
			leaves = []         
			if not node.is_leaf:             
				leaves += _get_leaves(node.left_child)             
				leaves += _get_leaves(node.right_child)         
			else:             
				leaves.append(node)         
			return leaves     
		return _get_leaves(root) 
			
	def _get_tree_rules(root):     
		def _get_rules(node):         
			rules = []         
			if node.rule is not None:             
				rules.append(node.rule)             
				rules += _get_rules(node.left_child)             
				rules += _get_rules(node.right_child)         
			return rules     
		return _get_rules(root)

class ProbabilisticTreeNode(TreeNode):     
	def predict(self, X):         
		"""         
		Multi-class predictions using the current node's rule         
		"""
		# Get probabilistic predictions         
		class_probabilities = self.predict_proba(X)
		
		# Find class with the max probability
		predictions = np.argmax(class_probabilities, axis=0)
		return predictions
		
		
	def predict_proba(self, X):
		"""         
		Probabilistic class predictions using the current node's rule         
		"""         
		class_probabilities = np.zeros((len(self.class_examples_idx.keys()), X.shape[0]))
			  
		# Push each example down the tree (an example is a row of X)         
		for i, x in enumerate(X):             
			x = x.reshape(1, -1)             
			current_node = self
					  
			# While we are not at a leaf             
			while not current_node.is_leaf:                 
			 # If the rule of the current node returns TRUE, branch left                 
			 if current_node.rule.classify(x):                     
				 current_node = current_node.left_child                 
			 # Otherwise, branch right                 
			 else:                     
				 current_node = current_node.right_child    
						  
			# A leaf has been reached. Use the leaf class proportions as the the class probabilities.
			for c in self.class_examples_idx.keys():
			  class_probabilities[c][i] = current_node.breiman_info.p_j_given_t[c]     
		                 
		return class_probabilities
