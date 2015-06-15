__author__ = 'Alexandre'

import h5py as h
import logging
import numpy as np
from sys import argv

def benjamini_yekutieli_correction(pvals, alpha=0.05):
    pvals = np.asarray(pvals)
    n_pvals = len(pvals)

    pvals_sortind = np.argsort(pvals)
    pvals_sorted = pvals[pvals_sortind]
    sortrevind = pvals_sortind.argsort()
    del pvals_sortind

    cm = np.sum(1./np.arange(1, n_pvals+1))

    # Find the greatest k such that the condition below is met. Reject null hypothesis for tests 1, ... k
    reject = pvals_sorted <= alpha * np.arange(1, n_pvals+1) / float(n_pvals) / cm
    if reject.any():
        rejectmax = max(np.nonzero(reject)[0])
    else:
        rejectmax = 0
    reject[:rejectmax] = True

    return reject[sortrevind]

def _auto_adjusting_alpha(mt_correction_func, p_values, alpha, max_significant):
    """
    Multiple testing correction with a maximum number of features that can be returned.

    Parameters:
    -----------
    mt_correction_func: function with arguments (p_values: ndarray of type float, alpha: float) that returns a boolean
                        mask (ndarray of the same shape as p_values)
        The multiple testing correction function.
    p_values: ndarray, dtype=float, shape=(n_p_values,)
        The p-values
    alpha: float
        The alpha parameter for the multiple testing correction function.
    max_significant: uint
        The maximum number of features that can be returned. If surpassed at a given value of alpha, alpha will be
        reduced until the number of significant features is less than or equal to max_significant.

    Returns:
    --------
    significant: ndarray, dtype=np.bool, shape=(n_p_values,)
        A boolean mask with the value TRUE for the features that are significant after correcting for multiple testing.
    alpha_final:
        The value of alpha that was used in order to satisfy computational requirements.
    """
    min_alpha = 1e-14
    max_alpha = 1000000.0
    alpha_decrease_factor = 2
    alpha_increase_factor = 2
    previous_action = None
    previous_alpha = None

    while min_alpha < alpha < max_alpha:
        significance_mask = mt_correction_func(p_values, alpha)
        n_significant = significance_mask.sum()

        if 0 < n_significant <= max_significant:
            break

        elif n_significant > max_significant:
            if previous_action is None or previous_action == "decrease":
                logging.debug("Too many significant features (%d) were found at alpha=%.14f. Decreasing alpha." % (n_significant, alpha))
                previous_alpha = alpha
                alpha /= alpha_decrease_factor
                previous_action = "decrease"
            else:
                # binary search mode
                significance_mask, alpha = _auto_adjusting_alpha_binary_search(mt_correction_func, p_values,
                                                                               previous_alpha, alpha, max_significant)
                break

        else:                                   # No significant features were found
            if previous_action is None or previous_action == "increase":
                logging.debug("No significant features (%d) were found at alpha=%.14f. Increasing alpha." % (n_significant, alpha))
                previous_alpha = alpha
                alpha *= alpha_increase_factor
                previous_action = "increase"
            else:
                # binary serach mode
                significance_mask, alpha = _auto_adjusting_alpha_binary_search(mt_correction_func, p_values, alpha,
                                                                               previous_alpha, max_significant)
                break

    logging.debug("%d significant features were found at alpha=%.14f" % (n_significant, alpha))
    return significance_mask, alpha

def _auto_adjusting_alpha_binary_search(mt_correction_func, p_values, alpha_min, alpha_max, max_significant):
    max_iter = 10
    bin_max = alpha_max
    bin_min = alpha_min

    n_iter = 0
    while n_iter < max_iter:
        mid = bin_min + (bin_max - bin_min) / 2

        significance_mask = mt_correction_func(p_values, mid)
        n_significant = significance_mask.sum()

        if 0 < n_significant <= max_significant:
            break
        elif 0 == n_significant:
            bin_min = mid
            logging.debug("No significant features (%d) were found at alpha=%.14f. Increasing alpha." % (n_significant,
                                                                                                         mid))
        else:
            bin_max = mid
            logging.debug("Too many significant features (%d) were found at alpha=%.14f. Decreasing alpha." % (n_significant,
                                                                                                               mid))

        n_iter += 1

    if n_significant == 0 or n_significant > max_significant:
        raise RuntimeError("Could find a proper alpha.")

    return significance_mask, mid


if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s.%(msecs)d %(levelname)s - [%(process)d] - %(module)s - %(funcName)s: %(message)s")

    # Arguments
    dataset_file = argv[1]
    alpha = float(argv[2])
    max_significant = int(argv[3])

    dataset = h.File(dataset_file)

    # Find all the groups in the file
    groups = [dataset["train"]] + [v for k, v in dataset["train/cv"].iteritems() if "fold_" in k]
    print "Found", len(groups), "groups."

    # Find the significant attributes in each group (presence attributes only and also do not consider attributes
    # that are in all examples)
    for group in groups:
        logging.debug("Processing group %s:" % str(group))

        n_attributes = group["p_values"].shape[0]
        type_by_attribute = dataset["type_by_attribute"][...]
        presence_attribute_selector = type_by_attribute == 0
        p_values = group["p_values"][presence_attribute_selector]

        presence_attribute_significance_mask, alpha_used = _auto_adjusting_alpha(benjamini_yekutieli_correction, p_values, alpha, max_significant)
        significance_mask = np.zeros(n_attributes, dtype=np.bool)
        significance_mask[presence_attribute_selector] = presence_attribute_significance_mask

        logging.debug("Saving the significance masks in datasets")
        mask_dataset = group.create_dataset("benjamini_yekutieli_significance_mask", data=significance_mask)
        mask_dataset.attrs["alpha_initial"] = alpha
        mask_dataset.attrs["alpha_used"] = alpha_used
        mask_dataset.attrs["max_significant"] = max_significant

    dataset.close()