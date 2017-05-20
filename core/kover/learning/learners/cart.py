#!/usr/bin/env python
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
from collections import defaultdict, deque 
from math import ceil 

from ..common.models import cart, CART_Model
from ..common.tree import ProbabilisticTreeNode

UTIL_BLOCK_SIZE = 1000000

class DecisionTreeClassifier(object):
	def __init__(self, criterion, max_depth, min_samples_split, class_importance):
		# Validate the node splitting criterion
		supported_criteria = ["gini"]
		if criterion not in supported_criteria:
			raise ValueError("The supporting splitting criteria are: %s." % str(supported_criteria))
		self.criterion = criterion	            

		# Validate the maximum depth of the tree
		if max_depth < 1:
			raise ValueError("The maximum tree depth must be greater than 1.")
		self.max_depth = max_depth
		
		# Validate the minimum number of examples required to split a node         
		# (This is a proportion of the number of training examples)
		
		if min_samples_split < 0.0 or min_samples_split > 1.0:
			raise ValueError( "The minimum number of examples (proportion) used to split a node must be between 0 and 1.")
		self.min_samples_split = float(min_samples_split)
		
		self.class_importance = class_importance
		
		self.rule_importances = defaultdict(float)
		self.model = CART_Model()
		
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
		n_total_class_examples = {c:float(len(ids)) for c, ids in examples_idx.items()}
		n_total_examples = sum(n_total_class_examples.values())
		
		# Compute the class priors based on the importance of making errors on each class.
		# See Section "4.4 Priors and Variable Misclassification Costs" in Breiman et al. 1984.
		priors = [1.0 * n_examples/ n_total_examples for n_examples in n_total_class_examples.values()]
		
		denom = sum([importance * prior for importance, prior in zip(self.class_importance, priors)])
		
		altered_priors = {c: 1.0 * importance * prior /denom for c, importance, prior in \
									zip(n_total_class_examples.keys(), self.class_importance, priors)}
		del n_total_examples, priors, denum
		
		# Criteria for node impurity and splitting
		def _gini_impurity(n_class_examples, multiply_by_node_proba=False):
			p_class_node = {c: 1.0*self.class_importance[c]*n_class_examples[c]/n_total_class_examples[c] \
											for c in n_class_examples.keys()}
			
			node_resubstitution_estimate = sum(p_class_node.values())
			p_class_given_node = {c: p_class_node[c]/ node_resubstitution_estimate for c in p_class_node.keys()}
			del p_class_node
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
						for c in example_idx.keys()}
			right_count = {c:np.asarray(example_idx[c].shape[0] - left_count[c], dtype=np.float) for c in left_count.keys()}
			
			# Left child:
			gini = _gini_impurity(left_count, multiply_by_node_proba=True) / node.breiman_info.p_t
			# Right child:
			gini += _gini_impurity(right_count, multiply_by_node_proba=True) / node.breiman_info.p_t
			
			# Don't consider rules that lead to empty nodes
			gini[sum(left_count.values()) == 0] = np.infty
			gini[sum(right_count.values())  == 0] = np.infty
			
			return gini
		
		if self.criterion == "gini":             
			get_criterion = _gini_impurity             
			score_rules = _gini_rule_score             
			choice_func = min             
			node_type = ProbabilisticTreeNode
			
		def _find_best_split(node):
			"""
			Selects the best split according to the splitting criterion
			
			"""
			example_idx = node.class_examples_idx
			
			# Score all the rules according to the criterion
			rules_criterion = score_rules(example_idx, node)
			
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
			
			return selected_rule_idx, best_rules_idx, left_child_example_idx_by_class, \                    
						right_child_example_idx_by_class
		
		logging.debug("Training start.")
		
		root = node_type(class_examples_idx=example_idx, depth=0,                          
			criterion_value=get_criterion(n_total_class_examples),                          
			class_priors=altered_priors,                          
			total_n_examples_by_class=n_total_class_examples)
		
		nodes_to_split = deque([root])
		runtime_infos = {}
		current_depth = -1
		min_samples_split = int(ceil(self.min_samples_split * sum(n_total_class_examples.values())))
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
            
            node.left_child = node_type(parent=node,
                                        class_examples_idx=left_child_example_idx_by_class, depth=node.depth + 1,
                                        criterion_value=get_criterion(left_child_example_idx_by_class),
                                        class_priors=altered_priors,
                                        total_n_examples_by_class=n_total_class_examples)
                                        
			node.right_child = node_type(parent=node,
                                        class_examples_idx=right_child_example_idx_by_class, depth=node.depth + 1,
                                        criterion_value=get_criterion(right_child_example_idx_by_class),
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
		
		# Save the model
        self.model.add_decision_tree(root)
        
        # Normalize the variable importances
        logging.debug("Normalizing the variable importances.")
        variable_importance_sum = sum(v for v in self.rule_importances.itervalues())
        self.rule_importances = dict((r, i / variable_importance_sum) for r, i in self.rule_importances.iteritems())

        logging.debug("Training finished.")
        
	def predict(self, X):
        if not self._is_fitted():
            raise RuntimeError("The classifier must be fitted before predicting.")
        return self.model.predict(X)

    def predict_proba(self, X):
        if not self._is_fitted():
            raise RuntimeError("The classifier must be fitted before predicting.")
        return self.model.predict_proba(X)

    def _is_fitted(self):
        return self.model is not None
