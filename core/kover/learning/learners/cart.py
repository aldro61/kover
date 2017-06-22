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
from copy import deepcopy

from ..common.models import cart, CARTModel
from ..common.tree import ProbabilisticTreeNode

UTIL_BLOCK_SIZE = 1000000

class DecisionTreeClassifier(object):
	def __init__(self, criterion, max_depth, min_samples_split, class_importance):
		# Validate the node splitting criterion
		supported_criteria = ["gini", "cross-entropy"]
		if criterion not in supported_criteria:
			raise ValueError("The supporting splitting criteria are: %s." % str(supported_criteria))
		self.criterion = criterion	            

		# Validate the maximum depth of the tree
		if max_depth < 1:
			raise ValueError("The maximum tree depth must be greater than 1.")
		self.max_depth = max_depth
		
		# Validate the minimum number of examples required to split a node         
		if min_samples_split < 2.0:
			raise ValueError( "The minimum number of examples used to split a node must be 2 or greater.")
		self.min_samples_split = int(min_samples_split)
		
		self.class_importance = class_importance
		self.rule_importances = defaultdict(float)
		
	def fit(self, rules, rule_classifications, example_idx, rule_blacklist=None,             
			tiebreaker=None, level_callback=None, split_callback=None):
				
		"""         
		Fits the decision tree classifier         
		"""
		
		if level_callback is None:
			level_callback = lambda x: None
		
		if split_callback is None:             
			split_callback = lambda x, y: None   
			      
		if tiebreaker is None:             
			tiebreaker = lambda x: x 
			
		# Store important information about the training set
		n_total_class_examples = {c:float(len(ids)) for c, ids in example_idx.items()}
		n_total_examples = sum(n_total_class_examples.values())
		
		# Compute the class priors based on the importance of making errors on each class.
		# See Section "4.4 Priors and Variable Misclassification Costs" in Breiman et al. 1984.
		priors = [1.0 * n_examples/ n_total_examples for n_examples in n_total_class_examples.values()]

		denum = sum([importance * prior for importance, prior in zip(self.class_importance.values(), priors)])
		altered_priors = {c: 1.0 * importance * prior / denum for c, importance, prior in \
									zip(n_total_class_examples.keys(), self.class_importance.values(), priors)}
									
		del n_total_examples, priors
		
		# Criteria for node impurity and splitting
		
		#############################################
		#				GINI IMPURITY				#		
		#############################################
		
		def _gini_impurity(n_class_examples, multiply_by_node_proba=False):
			p_class_node = {c: 1.0*altered_priors[c]*n_class_examples[c]/n_total_class_examples[c] \
											for c in n_class_examples.keys()}
			
			node_resubstitution_estimate = sum(p_class_node.values())
			with np.errstate(divide='ignore', invalid='ignore'):
				p_class_given_node = {c: np.divide(p_class_node[c], node_resubstitution_estimate) for c in p_class_node.keys()}

			diversity_index = sum([p_class_given_node[i]*p_class_given_node[j] \
									for i in p_class_given_node.keys() for j in p_class_given_node.keys() if i != j])
			
			return diversity_index * (node_resubstitution_estimate if multiply_by_node_proba else 1.0)
				
		def _gini_rule_score(example_idx, node):
			logging.debug("Scoring rules with the gini impurity strategy")             
			# A k-mer rule splits the examples in two groups: those that don't have the k-mer in their             
			# genome (left child) and those that have the k-mer in their genome (right child).             
			# XXX: We keep only the first half of the rule list and classifications, since sum_rows returns the counts             
			# XXX: for presence and absence rules, which we don't need here.             
			# TODO: We could have a parameter called return_absence=True, which avoid using twice the RAM.
			
			last_presence_rule_idx = int(1.0 * len(rules) / 2)  # Just because sum_rows returns pres/abs rules

			left_count = {c:np.asarray(rule_classifications.sum_rows(example_idx[c])[: last_presence_rule_idx],dtype=np.float) \
						for c in example_idx.keys() if example_idx[c].size}
			right_count = {c:np.asarray(example_idx[c].shape[0] - left_count[c], dtype=np.float) for c in left_count.keys()}

			# Left child:
			gini = _gini_impurity(left_count, multiply_by_node_proba=True) / node.breiman_info.p_t
			# Right child:
			gini += _gini_impurity(right_count, multiply_by_node_proba=True) / node.breiman_info.p_t
			
			# Don't consider rules that lead to empty nodes
			gini[sum(left_count.values()) == 0] = np.infty
			gini[sum(right_count.values())  == 0] = np.infty

			return gini
			
		#############################################
		#				CROSS-ENTROPY				#		
		#############################################
		
		def _cross_entropy(n_class_examples, multiply_by_node_proba=False):
			p_class_node = {c: 1.0*altered_priors[c]*n_class_examples[c]/n_total_class_examples[c] \
											for c in n_class_examples.keys()}

			node_resubstitution_estimate = sum(p_class_node.values())
			with np.errstate(divide='ignore', invalid='ignore'):
				p_class_given_node = {c: np.divide(p_class_node[c], node_resubstitution_estimate) for c in p_class_node.keys()}
				diversity_index = -1.0 * sum([np.nan_to_num(p_class_given_node[c]*np.log(p_class_given_node[c])) \
																						for c in p_class_given_node.keys()])
			return diversity_index * (node_resubstitution_estimate if multiply_by_node_proba else 1.0)
		
		def _cross_entropy_rule_score(example_idx, node):
			logging.debug("Scoring rules with the gini impurity strategy")             
			# A k-mer rule splits the examples in two groups: those that don't have the k-mer in their             
			# genome (left child) and those that have the k-mer in their genome (right child).             
			# XXX: We keep only the first half of the rule list and classifications, since sum_rows returns the counts             
			# XXX: for presence and absence rules, which we don't need here.             
			# TODO: We could have a parameter called return_absence=True, which avoid using twice the RAM.
			
			last_presence_rule_idx = int(1.0 * len(rules) / 2)  # Just because sum_rows returns pres/abs rules

			left_count = {c:np.asarray(rule_classifications.sum_rows(example_idx[c])[: last_presence_rule_idx],dtype=np.float) \
						for c in example_idx.keys() if example_idx[c].size}
			right_count = {c:np.asarray(example_idx[c].shape[0] - left_count[c], dtype=np.float) for c in left_count.keys()}

			# Left child:
			cross_entropy = _cross_entropy(left_count, multiply_by_node_proba=True) / node.breiman_info.p_t
			# Right child:
			cross_entropy += _cross_entropy(right_count, multiply_by_node_proba=True) / node.breiman_info.p_t
			
			# Don't consider rules that lead to empty nodes
			cross_entropy[sum(left_count.values()) == 0] = np.infty
			cross_entropy[sum(right_count.values())  == 0] = np.infty

			return cross_entropy
		
		if self.criterion == "gini":             
			get_criterion = _gini_impurity             
			score_rules = _gini_rule_score
			choice_func = min 
		elif self.criterion == "cross-entropy":
			get_criterion = _cross_entropy
			score_rules = _cross_entropy_rule_score
			choice_func = min             
		node_type = ProbabilisticTreeNode
			
		def _find_best_split(node):
			"""
			Selects the best split according to the splitting criterion
			
			"""
			example_idx = node.class_examples_idx
			
			# Score all the rules according to the criterion
			rules_criterion = score_rules(example_idx=example_idx, 
										  node=node)
			
			# Check if we find a split
			if (choice_func == min and min(rules_criterion) == np.infty) or (choice_func == max and max(rules_criterion) == -np.infty):
				return None, None, None, None
				
			# Tiebreaker to select a single rule in case of ties
			candidate_rules_idx = np.where(rules_criterion == choice_func(rules_criterion))[0]
			logging.debug("There are %d equivalent rules before the tiebreaker." % len(candidate_rules_idx))
			best_rules_idx = tiebreaker(candidate_rules_idx)
			logging.debug("The are %d equivalent rules after the tiebreaker." % len(best_rules_idx))
			selected_rule_idx = best_rules_idx[0]
			
			# Dispatch the examples to the children nodes based on the splitting rule's predictions
			rule_preds = rule_classifications.get_columns(selected_rule_idx)
			
			left_child_example_idx_by_class = {c:example_idx[c][rule_preds[example_idx[c]] == 1] for c in example_idx.keys()}
			right_child_example_idx_by_class = {c:example_idx[c][rule_preds[example_idx[c]] == 0] for c in example_idx.keys()}
			
			return selected_rule_idx, best_rules_idx, left_child_example_idx_by_class, right_child_example_idx_by_class
		
		logging.debug("Training start.")
		
		root = node_type(class_examples_idx=example_idx, 
						 depth=0,                          
						 criterion_value=get_criterion(n_total_class_examples),                          
						 class_priors=altered_priors,                          
						 total_n_examples_by_class=n_total_class_examples)

		nodes_to_split = deque([root])
		runtime_infos = {}
		current_depth = -1
		min_samples_split = self.min_samples_split
		if min_samples_split < 2:
			min_samples_split = 2
		
		while len(nodes_to_split) > 0:
			node = nodes_to_split.popleft()
			
			# Check if we have reached a new depth
			if node.depth != current_depth:
				current_depth = node.depth

				runtime_infos["depth"] = current_depth

				logging.debug("The tree depth is %d" % current_depth)
				if current_depth > 0:
					# The level callback is called when all the nodes of a level have been created
					level_callback(runtime_infos)
				if current_depth == self.max_depth:
					logging.debug("The maximum tree depth has been reached. No more leaves will be split.")
					break  # We have reached the nodes of the last level, which must remain leaves
			
			# Check if the node to split is a pure leaf
			if 1.0 in node.class_proportions.values():
				logging.debug("The leaf is pure. It will not be split.")
				continue
				
			# Check if the HP constraints allows us to split this node
			if node.n_examples < min_samples_split:
				logging.debug("The leaf contains less examples (%d) than the minimum required to split (%d) a node. "                               
				"It will not be split." % (node.n_examples, min_samples_split))
				continue
				
			# Find the best rule to split the node
			selected_rule_idx, \
			equivalent_rule_idx, \
			left_child_example_idx_by_class, \
			right_child_example_idx_by_class = _find_best_split(node)
			
			# If we were incapable of splitting the node into two non-empty leafs
			if selected_rule_idx is None:
				logging.debug("Found no rule to split the node. The node will not be split.")
				continue
				
			# Perform the split
			node.rule = rules[selected_rule_idx]
			left_child_n_class_example = {c:len(idx) for c, idx in left_child_example_idx_by_class.items()}
			right_child_n_class_example = {c:len(idx) for c, idx in right_child_example_idx_by_class.items()}
			
			node.left_child = node_type(parent=node,
										class_examples_idx=left_child_example_idx_by_class, 
										depth=node.depth + 1,
										criterion_value=get_criterion(left_child_n_class_example),
										class_priors=altered_priors,
										total_n_examples_by_class=n_total_class_examples)
										
			node.right_child = node_type(parent=node,
										class_examples_idx=right_child_example_idx_by_class, 
										depth=node.depth + 1,
										criterion_value=get_criterion(right_child_n_class_example),
										class_priors=altered_priors,
										total_n_examples_by_class=n_total_class_examples)
										
			# Update rule importances
			# TODO: account for node probabilities (this is currently incorrect)
			self.rule_importances[str(node.rule)] += node.criterion_value - \
													 node.left_child.criterion_value - \
													 node.right_child.criterion_value
													 
			logging.debug("Split with rule %s." % node.rule)
			split_callback(node, equivalent_rule_idx)

			# Add the new child nodes to the splitting queue
			nodes_to_split.append(node.left_child)
			nodes_to_split.append(node.right_child)

			# Update the model in the runtime informations
			runtime_infos["model"] = root
			
		logging.debug("Done building the tree.")
		
		# Save the decision tree
		self.decision_tree = root
		
		# Normalize the variable importances
		logging.debug("Normalizing the variable importances.")
		variable_importance_sum = sum(v for v in self.rule_importances.itervalues())
		self.rule_importances = dict((r, i / variable_importance_sum) for r, i in self.rule_importances.iteritems())

		logging.debug("Training finished.")
		
	def predict(self, X):
		if not self._is_fitted():
			raise RuntimeError("The classifier must be fitted before predicting.")
		return self.decision_tree.predict(X)

	def predict_proba(self, X):
		if not self._is_fitted():
			raise RuntimeError("The classifier must be fitted before predicting.")
		return self.decision_tree.predict_proba(X)

	def _is_fitted(self):
		return self.decision_tree is not None
		
def _prune_tree(tree):
	"""
	Prunes a decision tree
	
	"""
	def _initial_pruning(root):
		# Get all nodes that have two leaves as children
		def _get_leaf_parents(root):
			def __get_leaf_parents(node):
				leaf_parents = []
				if not node.is_leaf:
					if node.left_child.is_leaf and node.right_child.is_leaf:
						leaf_parents.append(node)
					else:
						leaf_parents += __get_leaf_parents(node.left_child)
						leaf_parents += __get_leaf_parents(node.right_child)
				return leaf_parents
			return __get_leaf_parents(root)
			
		# Perform the initial pruning (Tmax -> T1)
		def __initial_pruning(parents):
			if len(parents) == 0:
				return
			node = parents.pop()
			if np.allclose(node.breiman_info.R_t, node.left_child.breiman_info.R_t + node.right_child.breiman_info.R_t):
				logging.debug("Converting node %s into a leaf" % node)
				del node.rule
				node.rule = None
				del node.left_child
				node.left_child = None
				del node.right_child
				node.right_child = None
				if not node.is_root and node.parent.left_child.is_leaf and node.parent.right_child.is_leaf:
					logging.debug("Adding the new leaf's parent to the list of leaf parents")
					parents.append(node.parent)
			__initial_pruning(parents)
			
		# Start the initial pruning
		__initial_pruning(_get_leaf_parents(root))
		
	def _find_weakest_links(root):
		def __find_weakest_links(node):
			if node.is_leaf:
				return np.inf, [node]
			else:
				RTt = sum(l.breiman_info.R_t for l in node.leaves)
				current_gt = float(node.breiman_info.R_t - RTt) / (len(node.leaves) - 1)
				left_min_gt, left_weakest_links = __find_weakest_links(node.left_child)
				right_min_gt, right_weakest_links = __find_weakest_links(node.right_child)
				
				if np.allclose(current_gt, min(left_min_gt, right_min_gt)):
					if np.allclose(left_min_gt, right_min_gt):
						return current_gt, [node] + left_weakest_links + right_weakest_links
					else:
						return current_gt, [node] + (left_weakest_links if left_min_gt < right_min_gt else right_weakest_links)
				elif current_gt < min(left_min_gt, right_min_gt):                         
					return current_gt, [node]                     
				elif np.allclose(left_min_gt, right_min_gt):                         
					return left_min_gt, left_weakest_links + right_weakest_links                     
				elif left_min_gt > right_min_gt:                         
					return right_min_gt, right_weakest_links                     
				elif left_min_gt < right_min_gt:                         
					return left_min_gt, left_weakest_links                     
				else:                         
					raise Exception("Unhandled case detected!")
						
		return __find_weakest_links(root)
			
	def _sequential_prune(root):
		Tmax = deepcopy(root)
		logging.debug("Initial pruning (Tmax >> T1)")
		_initial_pruning(Tmax)
		T1 = Tmax
		del Tmax
		
		def __sequential_prune(root):
			root = deepcopy(root)

			# Find the weakest link in the tree                 
			logging.debug("Computing link strenghts for each node")             
			min_gt, weakest_links = _find_weakest_links(root)
			
			# Prune the weakest link (and save alpha)
			logging.debug("Pruning occurs at alpha %.9f" % min_gt)
			# TODO: faire une fonction qui descend jusqu'aux feuilles et prune en remontant car perte de mémoire ici
			
			for n in weakest_links:                     
				del n.rule                     
				n.rule = None                     
				del n.left_child                     
				n.left_child = None                     
				del n.right_child                     
				n.right_child = None
				
			# Repeat until only the root node remains
			return [(min_gt, root)] + (__sequential_prune(root) if not root.is_leaf else [])
		
		logging.debug("Pruning sequentially until only the root remains (T1 >> ... >> {root}")
		return [(0, T1)] + __sequential_prune(T1)
		
	logging.debug("Copying model")
	tree = deepcopy(tree) 
	
	logging.debug("Generating the sequence of pruned trees")
	alphas, trees = zip(*_sequential_prune(tree))
	
	return alphas, trees
