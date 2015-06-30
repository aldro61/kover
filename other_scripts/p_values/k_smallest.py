__author__ = 'Alexandre'

import h5py as h
import logging
import numpy as np
from sys import argv



if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s.%(msecs)d %(levelname)s - [%(process)d] - %(module)s - %(funcName)s: %(message)s")

    # Arguments
    dataset_file = argv[1]
    n_to_keep = int(argv[2])

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
        presence_attribute_idx = np.where(type_by_attribute == 0)[0]
        p_values = group["p_values"][presence_attribute_idx]

        logging.debug("Sorting the p-values for the presence attributes")
        argsort = np.argsort(p_values)

        logging.debug("Keeping the %d smallest p-values" % n_to_keep)
        significance_mask = np.zeros(n_attributes, dtype=np.bool)
        significance_mask[presence_attribute_idx[argsort[: n_to_keep]]] = True

        logging.debug("Saving the significance masks in datasets")
        mask_dataset = group.create_dataset("%d_smallest_significance_mask" % n_to_keep, data=significance_mask)

    dataset.close()