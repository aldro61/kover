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
import logging
import numpy as np


from .models import BaseModel


class BreimanInfo(object):
	def __init__(self, node_n_examples_by_class, class_priors, total_n_examples_by_class):

		# Eq. 2.2 Probability that an example is in class j and falls into node t
		self.p_j_t = {c: class_priors[c] * node_n_examples_by_class[c] / total_n_examples_by_class[c]
					     for c in class_priors.iterkeys()}

		# Eq. 2.3 Probability that any example falls in node t
		self.p_t = sum(self.p_j_t.values())

		# Eq. 2.4 Probability that an example is in class j given that it falls in node t
		self.p_j_given_t = {c: self.p_j_t[c] / self.p_t for c in class_priors.iterkeys()}

		# Def. 2.10 Probability of misclassification given that an example falls into node t
		self.r_t = 1.0 - max(self.p_j_given_t.values())

		# Contribution of the node to the tree's overall missclassification rate
		self.R_t = self.r_t * self.p_t


class TreeNode(BaseModel):
	# TODO: implement memoizing for attributes that require longer computation
	def __init__(self, depth, class_examples_idx, total_n_examples_by_class, class_priors, rule=None, parent=None,
				 left_child=None, right_child=None, criterion_value=-1):
		self.rule = rule
		self.parent = parent
		self.left_child = left_child
		self.right_child = right_child
		self.class_examples_idx = class_examples_idx
		self.depth = depth
		self.criterion_value = criterion_value

		assert isinstance(class_examples_idx, dict)
		assert isinstance(total_n_examples_by_class, dict)
		assert isinstance(class_priors, dict)

		n_examples_by_class = {c: len(c_idx) for c, c_idx in class_examples_idx.iteritems()}
		self.breiman_info = BreimanInfo(node_n_examples_by_class=n_examples_by_class,
										class_priors=class_priors,
										total_n_examples_by_class=total_n_examples_by_class)

	@property
	def is_leaf(self):
		"""
		Returns true if current node is a leaf

		"""
		return self.rule is None and self.left_child is None and self.right_child is None

	@property
	def is_root(self):
		"""
		Returns true if current node is the tree root

		"""
		return self.parent is None

	@property
	def n_examples(self):
		"""
		Returns the number of examples in the node

		"""
		return  sum(len(c_idx) for c_idx in self.class_examples_idx.itervalues())

	@property
	def class_proportions(self):
		"""
		Returns the proportion of examples of each class in the node

		"""
		n_examples = self.n_examples
		return {c: float(len(c_idx)) / n_examples for c, c_idx in self.class_examples_idx.iteritems()}

	@property
	def class_prediction(self):
		"""
		Returns the class predicted by this node as a leaf

		"""
		return self.breiman_info.p_j_given_t.keys()[np.argmax(self.breiman_info.p_j_given_t.values())]

	@property
	def rules(self):
		"""
		Returns all the rules of the tree

		"""
		return _get_tree_rules(self)

	@property
	def leaves(self):
		"""
		Returns all the leaves of the tree

		"""
		return _get_tree_leaves(self)

	@property
	def tree_depth(self):
		"""
		Returns the depth of the tree

		"""
		return _get_tree_depth(self)

	def __iter__(self):
		"""
		Yields all the rules of the tree in preorder

		"""
		def _preorder(node):
			nodes = [node]
			if not node.is_leaf:
				nodes += _preorder(node.left_child)
				nodes += _preorder(node.right_child)
			return nodes
		nodes = _preorder(self)

		for node_id, node in zip(range(len(nodes)), nodes):
			yield node_id, node

	def __len__(self):
		"""
		Returns the number of nodes in the tree

		"""
		return len(self.rules) + len(self.leaves)

	def __str__(self, depth=0):
		"""
		Convert the tree to a text representation

		"""
		tree_str = ""

		# Case : node is a leaf
		if self.is_leaf:
			tree_str += "\n" + ("    " * depth) + str(self.class_prediction)

		# Case : node has two children
		else:
			# Print right branch
			tree_str += self.right_child.__str__(depth=depth + 1)

			# Print own value
			tree_str += "\n" + ("    " * depth + "   ") + str("/")
			tree_str += "\n" + ("    " * depth) + str(self.rule)
			tree_str += "\n" + ("    " * depth + "   ") + str("\\")

			# Print left_child branch
			tree_str += self.left_child.__str__(depth=depth + 1)

		return tree_str


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
		if not node.is_leaf:
			rules.append(node.rule)
			rules += _get_rules(node.left_child)
			rules += _get_rules(node.right_child)
		return rules
	return _get_rules(root)


def _get_tree_depth(root):
	"""
	Assumes that the depth is builtin the tree nodes, does not calculate it.

	"""
	def _get_depth(node):
		if not node.is_leaf:
			depth = max(_get_depth(node.left_child), _get_depth(node.right_child))
		else:
			depth = node.depth
		return depth
	return _get_depth(root)


class ProbabilisticTreeNode(TreeNode):
	def predict(self, X):
		"""
		Multi-class predictions using the current node's rule

		"""
		X = np.ascontiguousarray(X)

		# Get probabilistic predictions
		class_probabilities = self.predict_proba(X)

		# Find class with the max probability
		predictions = np.argmax(class_probabilities, axis=0)
		return predictions


	def predict_proba(self, X):
		"""
		Probabilistic class predictions using the current node's rule

		"""
		X = np.ascontiguousarray(X)

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
			for c in self.class_examples_idx.iterkeys():
			  class_probabilities[c][i] = current_node.breiman_info.p_j_given_t[c]

		return class_probabilities
