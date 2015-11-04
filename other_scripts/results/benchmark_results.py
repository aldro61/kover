#!/usr/bin/env python
# -*- coding:utf-8 -*-

import cPickle as c
import fnmatch
import os
from sys import argv

import h5py as h
import numpy as np

def find_results(dir):
    results = []
    for root, dirnames, filenames in os.walk(dir):
        for filename in fnmatch.filter(filenames, 'dOut'):
            f = open(os.path.join(root, filename))
            try:
                while True:
                    results.append(c.load(f))
            except:
                pass
    return results

def get_immediate_subdirectories(a_dir):
    return [name for name in os.listdir(a_dir)
            if os.path.isdir(os.path.join(a_dir, name))]

def get_scm_results(results):
    print "Finding splits..."
    splits = sorted(list(set([r["config"]["dataset_file"] for r in results])))
    print "Found", len(splits), "splits."
    print

    best_cv_score_by_split = []
    best_train_risk_by_split = []
    best_test_risk_by_split = []
    best_test_f1score_by_split = []
    best_test_precision_by_split = []
    best_test_recall_by_split = []
    best_test_sensitivity_by_split = []
    best_test_specificity_by_split = []
    best_number_of_attributes_by_split = []

    for split in splits:
        split_res_counter = 0
        print
        print "Split:", split
        print "*" * 100

        best_result_cv_score = 1.0
        best_result = None
        best_result_config = None
        best_model_size = -1

        for result in results:

            # Only look at the results for this split
            if result["config"]["dataset_file"] != split:
                continue

            split_res_counter += 1

            config = result["config"]
            cv_score = result["cv"]["best"]["metrics_test"]["risk"][0]

            if cv_score < best_result_cv_score:
                best_result_cv_score = cv_score
                best_result = result
                best_result_config = config
                best_model_size = result["cv"]["best"]["average_model_complexity"]
            elif cv_score == best_result_cv_score:
                if result["cv"]["best"]["average_model_complexity"] < best_model_size:
                    best_result_cv_score = cv_score
                    best_result = result
                    best_result_config = config
                    best_model_size = result["cv"]["best"]["average_model_complexity"]
                elif result["cv"]["best"]["average_model_complexity"] == best_model_size and np.abs(1.0 - config["scm_p"]) < np.abs(1.0 - best_result_config["scm_p"]):
                    best_result_cv_score = cv_score
                    best_result = result
                    best_result_config = config
                    best_model_size = result["cv"]["best"]["average_model_complexity"]
                    
        print "The split had %d results" % split_res_counter
        print best_result_config
        print

        best_cv_score_by_split.append(best_result_cv_score)
        best_train_risk_by_split.append(best_result["full"]["metrics_train"]["risk"][0])
        best_test_risk_by_split.append(best_result["full"]["metrics_test"]["risk"][0])
        best_test_f1score_by_split.append(best_result["full"]["metrics_test"]["f1_score"][0])
        best_test_precision_by_split.append(best_result["full"]["metrics_test"]["precision"][0])
        best_test_recall_by_split.append(best_result["full"]["metrics_test"]["recall"][0])
        best_test_sensitivity_by_split.append(best_result["full"]["metrics_test"]["sensitivity"][0])
        best_test_specificity_by_split.append(best_result["full"]["metrics_test"]["specificity"][0])
        best_number_of_attributes_by_split.append(best_result["full"]["algorithm_specific"]["n_attributes"])

    return {"train_risk_mean": np.mean(best_train_risk_by_split),
            "train_risk_std": np.std(best_train_risk_by_split),
            "train_accuracy_mean": np.mean(1.0 - np.array(best_train_risk_by_split, dtype=np.float)),
            "train_accuracy_std": np.std(1.0 - np.array(best_train_risk_by_split, dtype=np.float)),
            "test_risk_mean": np.mean(best_test_risk_by_split),
            "test_risk_std": np.std(best_test_risk_by_split),
            "test_accuracy_mean": np.mean(1.0 - np.array(best_test_risk_by_split, dtype=np.float)),
            "test_accuracy_std": np.std(1.0 - np.array(best_test_risk_by_split, dtype=np.float)),
            "test_sensitivity_mean": np.mean(best_test_sensitivity_by_split),
            "test_sensitivity_std": np.std(best_test_sensitivity_by_split),
            "test_specificity_mean": np.mean(best_test_specificity_by_split),
            "test_specificity_std": np.std(best_test_specificity_by_split),
            "complexity_mean": np.mean(best_number_of_attributes_by_split),
            "complexity_std": np.std(best_number_of_attributes_by_split)}

def get_chi2_ksmall_decision_tree_results(results):
    print "Finding splits..."
    splits = sorted(list(set([r["config"]["dataset_file"] for r in results])))
    print "Found", len(splits), "splits."
    print

    best_cv_score_by_split = []
    best_train_risk_by_split = []
    best_test_risk_by_split = []
    best_test_f1score_by_split = []
    best_test_precision_by_split = []
    best_test_recall_by_split = []
    best_test_sensitivity_by_split = []
    best_test_specificity_by_split = []
    best_number_of_rules_by_split = []

    for split in splits:
        split_res_counter = 0
        print
        print "Split:", split
        print "*" * 100

        best_result_cv_score = 1.0
        best_result = None
        best_result_config = None
        best_depth = -1

        for result in results:

            # Only look at the results for this split
            if result["config"]["dataset_file"] != split:
                continue

            split_res_counter += 1

            config = result["config"]
            cv_score = result["cv"]["best"]["metrics_test"]["risk"][0]

            if cv_score < best_result_cv_score:
                best_result_cv_score = cv_score
                best_result = result
                best_result_config = config
                best_depth = config["decision_tree_depth"]
            elif cv_score == best_result_cv_score:
                if config["decision_tree_depth"] < best_depth:
                    best_result_cv_score = cv_score
                    best_result = result
                    best_result_config = config
                    best_depth = config["decision_tree_depth"]
        print "The split had %d results" % split_res_counter
        print best_result_config
        print

        best_cv_score_by_split.append(best_result_cv_score)
        best_train_risk_by_split.append(best_result["full"]["metrics_train"]["risk"][0])
        best_test_risk_by_split.append(best_result["full"]["metrics_test"]["risk"][0])
        best_test_f1score_by_split.append(best_result["full"]["metrics_test"]["f1_score"][0])
        best_test_precision_by_split.append(best_result["full"]["metrics_test"]["precision"][0])
        best_test_recall_by_split.append(best_result["full"]["metrics_test"]["recall"][0])
        best_test_sensitivity_by_split.append(best_result["full"]["metrics_test"]["sensitivity"][0])
        best_test_specificity_by_split.append(best_result["full"]["metrics_test"]["specificity"][0])
        best_number_of_rules_by_split.append(best_result["full"]["algorithm_specific"]["n_rules"])

    return {"train_risk_mean": np.mean(best_train_risk_by_split),
            "train_risk_std": np.std(best_train_risk_by_split),
            "train_accuracy_mean": np.mean(1.0 - np.array(best_train_risk_by_split, dtype=np.float)),
            "train_accuracy_std": np.std(1.0 - np.array(best_train_risk_by_split, dtype=np.float)),
            "test_risk_mean": np.mean(best_test_risk_by_split),
            "test_risk_std": np.std(best_test_risk_by_split),
            "test_accuracy_mean": np.mean(1.0 - np.array(best_test_risk_by_split, dtype=np.float)),
            "test_accuracy_std": np.std(1.0 - np.array(best_test_risk_by_split, dtype=np.float)),
            "test_sensitivity_mean": np.mean(best_test_sensitivity_by_split),
            "test_sensitivity_std": np.std(best_test_sensitivity_by_split),
            "test_specificity_mean": np.mean(best_test_specificity_by_split),
            "test_specificity_std": np.std(best_test_specificity_by_split),
            "complexity_mean": np.mean(best_number_of_rules_by_split),
            "complexity_std": np.std(best_number_of_rules_by_split)}

def get_chi2_ksmall_l1svm_results(results):
    print "Finding splits..."
    splits = sorted(list(set([r["config"]["dataset_file"] for r in results])))
    print "Found", len(splits), "splits."
    print

    best_cv_score_by_split = []
    best_train_risk_by_split = []
    best_test_risk_by_split = []
    best_test_f1score_by_split = []
    best_test_precision_by_split = []
    best_test_recall_by_split = []
    best_test_sensitivity_by_split = []
    best_test_specificity_by_split = []
    best_number_of_features_by_split = []

    for split in splits:
        split_res_counter = 0
        print
        print "Split:", split
        print "*" * 100

        best_result_cv_score = 1.0
        best_result = None
        best_result_config = None
        best_c = -1

        for result in results:

            # Only look at the results for this split
            if result["config"]["dataset_file"] != split:
                continue

            split_res_counter += 1

            config = result["config"]
            cv_score = result["cv"]["best"]["metrics_test"]["risk"][0]

            if cv_score < best_result_cv_score:
                best_result_cv_score = cv_score
                best_result = result
                best_result_config = config
                best_c = config["svm_c"]
            elif cv_score == best_result_cv_score:
                if config["svm_c"] < best_c:
                    best_result_cv_score = cv_score
                    best_result = result
                    best_result_config = config
                    best_c = config["svm_c"]
        print "The split had %d results" % split_res_counter
        print best_result_config
        print

        best_cv_score_by_split.append(best_result_cv_score)
        best_train_risk_by_split.append(best_result["full"]["metrics_train"]["risk"][0])
        best_test_risk_by_split.append(best_result["full"]["metrics_test"]["risk"][0])
        best_test_f1score_by_split.append(best_result["full"]["metrics_test"]["f1_score"][0])
        best_test_precision_by_split.append(best_result["full"]["metrics_test"]["precision"][0])
        best_test_recall_by_split.append(best_result["full"]["metrics_test"]["recall"][0])
        best_test_sensitivity_by_split.append(best_result["full"]["metrics_test"]["sensitivity"][0])
        best_test_specificity_by_split.append(best_result["full"]["metrics_test"]["specificity"][0])
        best_number_of_features_by_split.append(best_result["full"]["algorithm_specific"]["n_features"])

    return {"train_risk_mean": np.mean(best_train_risk_by_split),
            "train_risk_std": np.std(best_train_risk_by_split),
            "train_accuracy_mean": np.mean(1.0 - np.array(best_train_risk_by_split, dtype=np.float)),
            "train_accuracy_std": np.std(1.0 - np.array(best_train_risk_by_split, dtype=np.float)),
            "test_risk_mean": np.mean(best_test_risk_by_split),
            "test_risk_std": np.std(best_test_risk_by_split),
            "test_accuracy_mean": np.mean(1.0 - np.array(best_test_risk_by_split, dtype=np.float)),
            "test_accuracy_std": np.std(1.0 - np.array(best_test_risk_by_split, dtype=np.float)),
            "test_sensitivity_mean": np.mean(best_test_sensitivity_by_split),
            "test_sensitivity_std": np.std(best_test_sensitivity_by_split),
            "test_specificity_mean": np.mean(best_test_specificity_by_split),
            "test_specificity_std": np.std(best_test_specificity_by_split),
            "complexity_mean": np.mean(best_number_of_features_by_split),
            "complexity_std": np.std(best_number_of_features_by_split)}

def get_chi2_ksmall_l2svm_results(results):
    return get_chi2_ksmall_l1svm_results(results)

def get_chi2_ksmall_scm_results(results):
    return get_scm_results(results)

def get_baseline_results(results):
    print "Finding splits..."
    splits = sorted(list(set([r["config"]["dataset_file"] for r in results])))
    print "Found", len(splits), "splits."
    print

    best_cv_score_by_split = []
    best_train_risk_by_split = []
    best_test_risk_by_split = []
    best_test_f1score_by_split = []
    best_test_precision_by_split = []
    best_test_recall_by_split = []
    best_test_sensitivity_by_split = []
    best_test_specificity_by_split = []

    for split in splits:
        split_res_counter = 0
        print
        print "Split:", split
        print "*" * 100

        best_result_cv_score = 1.0
        best_result = None
        best_result_config = None
        best_c = -1

        for result in results:

            # Only look at the results for this split
            if result["config"]["dataset_file"] != split:
                continue

            split_res_counter += 1

            config = result["config"]
            cv_score = result["cv"]["best"]["metrics_test"]["risk"][0]

            if cv_score < best_result_cv_score:
                best_result_cv_score = cv_score
                best_result = result
                best_result_config = config

        print "The split had %d results" % split_res_counter
        print best_result_config
        print

        best_cv_score_by_split.append(best_result_cv_score)
        best_train_risk_by_split.append(best_result["full"]["metrics_train"]["risk"][0])
        best_test_risk_by_split.append(best_result["full"]["metrics_test"]["risk"][0])
        best_test_f1score_by_split.append(best_result["full"]["metrics_test"]["f1_score"][0])
        best_test_precision_by_split.append(best_result["full"]["metrics_test"]["precision"][0])
        best_test_recall_by_split.append(best_result["full"]["metrics_test"]["recall"][0])
        best_test_sensitivity_by_split.append(best_result["full"]["metrics_test"]["sensitivity"][0])
        best_test_specificity_by_split.append(best_result["full"]["metrics_test"]["specificity"][0])

    return {"train_risk_mean": np.mean(best_train_risk_by_split),
            "train_risk_std": np.std(best_train_risk_by_split),
            "train_accuracy_mean": np.mean(1.0 - np.array(best_train_risk_by_split, dtype=np.float)),
            "train_accuracy_std": np.std(1.0 - np.array(best_train_risk_by_split, dtype=np.float)),
            "test_risk_mean": np.mean(best_test_risk_by_split),
            "test_risk_std": np.std(best_test_risk_by_split),
            "test_accuracy_mean": np.mean(1.0 - np.array(best_test_risk_by_split, dtype=np.float)),
            "test_accuracy_std": np.std(1.0 - np.array(best_test_risk_by_split, dtype=np.float)),
            "test_sensitivity_mean": np.mean(best_test_sensitivity_by_split),
            "test_sensitivity_std": np.std(best_test_sensitivity_by_split),
            "test_specificity_mean": np.mean(best_test_specificity_by_split),
            "test_specificity_std": np.std(best_test_specificity_by_split),
            "complexity_mean": 0.0,
            "complexity_std": 0.0}
   

if __name__ == "__main__":
    species = ["cdiff", "mtuberculosis", "paeruginosa", "spneumoniae"]
    specie_printable_names = ["C. difficile", "M. tuberculosis", "P. aeruginosa", "S. pneumoniae"]
    our_method = "scm"
    methods = ["scm", "chi2_ksmall_decision_tree", "chi2_ksmall_l1svm", "chi2_ksmall_l2svm", "baseline"]
    method_printable_names = ["SCM", "$\\chi^2 + \\mbox{CART}$", "$\\chi^2 + \\mbox{L1SVM}$", "$\\chi^2 + \\mbox{L2SVM}$", "Baseline"]
    metric_1 = "test_accuracy"
    metric_1_decimals = 3
    metric_2 = "complexity"
    metric_2_decimals = 1
    include_second_metric = True
    include_standard_deviation = False  # The naming convention is metric_mean and metric_std the code will try to get these.

    table = "\\begin{tabular}{l|" + "l" * len(methods) + "} \n"
    table += "\\hline \n"
    table += "Dataset & " + " & ".join([m.replace("_", "\\_") for m in method_printable_names]) + "\\\\ \n"

    from collections import defaultdict
    method_results = defaultdict(lambda: defaultdict(dict))
    for specie, specie_printable_name in zip(species, specie_printable_names):
        if os.path.exists(specie):
            datasets = sorted(get_immediate_subdirectories(specie))
            print "Specie:", specie
            table += "\\hline \n"
            table += "\\textbf{\\textit{%s}}" % specie_printable_name + " & " * len(methods) + "\\\\ \n"
            for dataset in datasets:
                print "Dataset:", dataset

                dataset_table_line = "\\,\\, " + dataset.title() + " & "

                method_printable_results = []
                for method in methods:
                    print "Method:", method
                    print "Searching for results in:", "./%s/%s/benchmark/%s"%(specie, dataset, method)
                    results = []
                    if os.path.exists("./%s/%s/benchmark/%s"%(specie, dataset, method)):
                        results = find_results("./%s/%s/benchmark/%s"%(specie, dataset, method))

                    if len(results) == 0:
                        print "Warning: No results found for %s %s %s" % (specie, dataset, method)
                        printable_result = "NA (NA)"
                    else:
                        # Find the best result
                        results = globals()["get_%s_results" % method](results)
                        if include_second_metric:
                            if include_standard_deviation:
                            	printable_result = "%%.%df; %%.%df (%%.%df; %%.%df)" % (metric_1_decimals, metric_1_decimals, metric_2_decimals, metric_2_decimals) % (round(results[metric_1 + "_mean"], metric_1_decimals), round(results[metric_1 + "_std"], metric_1_decimals), round(results[metric_2 + "_mean"], metric_2_decimals), round(results[metric_2 + "_std"], metric_2_decimals))
                            else:
                            	printable_result = "%%.%df (%%.%df)" % (metric_1_decimals, metric_2_decimals)  % (round(results[metric_1 + "_mean"], metric_1_decimals), round(results[metric_2 + "_mean"], metric_2_decimals))
                        else:
                            if include_standard_deviation:
                                printable_result = "%%.%df; %%.%df" % (metric_1_decimals, metric_1_decimals) % (round(results[metric_1 + "_mean"], metric_1_decimals), round(results[metric_1 + "_std"], metric_1_decimals))
                            else:
                                printable_result = "%%.%df" % (metric_1_decimals)  % (round(results[metric_1 + "_mean"], metric_1_decimals))
                        method_results[method][specie + "_" + dataset] = results

                    method_printable_results.append(printable_result)

                dataset_table_line += " & ".join(method_printable_results) + "\\\\ \n"
                table += dataset_table_line
        else:
            print "Warning: No experiments found for specie %s"%specie

    # Get the p-value for the comparison to our method for each metric
    p_values_metric_1 = []
    p_values_metric_2 = []
    our_method_results = set(method_results[our_method].keys())
    for m in methods:
        if m == our_method:
            p_values_metric_1.append("")
            p_values_metric_2.append("")
        else:
            # Find all the datasets for which there are results between m and our method
            common_results = our_method_results.intersection(set(method_results[m].keys()))

            if len(common_results) != len(our_method_results):
                print "WARNING: Could not compute p-values for methods %s and %s since they do not have results for the same datasets" % (our_method, m)
                p_values_metric_1.append("NA")
                p_values_metric_2.append("NA")
            else:
                # Compute the p-values for both metrics
                our_results_metric_1 = [method_results[our_method][dataset][metric_1 + "_mean"] for dataset in common_results]
                our_results_metric_2 = [method_results[our_method][dataset][metric_2 + "_mean"] for dataset in common_results]
                other_results_metric_1 = [method_results[m][dataset][metric_1 + "_mean"] for dataset in common_results]
                other_results_metric_2 = [method_results[m][dataset][metric_2 + "_mean"] for dataset in common_results]

                from scipy.stats import wilcoxon
                p_value_metric_1 = round(wilcoxon(our_results_metric_1, other_results_metric_1, zero_method="wilcox")[1], 4)
                p_value_metric_2 = round(wilcoxon(our_results_metric_2, other_results_metric_2, zero_method="wilcox")[1], 4)

                p_values_metric_1.append(p_value_metric_1)
                p_values_metric_2.append(p_value_metric_2)

    table += metric_1.replace("_", "\\_") + " p-values: & " + " & ".join([str(v) for v in p_values_metric_1]) + "\\\\ \n"
    table += metric_2.replace("_", "\\_") + " p-values: & " + " & ".join([str(v) for v in p_values_metric_2]) + "\\\\ \n"

    table += "\\end{tabular}"
    print table

    # scm_results = get_scm_results()
    # chi2_kbest_dt_results = get_chi2_kbest_decision_tree_results()
    # chi2_kbest_l1svm_results = get_chi2_kbest_l1svm_results()
    # chi2_kbest_l1svm_results = get_chi2_kbest_l2svm_results()
