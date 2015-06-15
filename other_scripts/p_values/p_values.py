"""
Computes the p-values for a dataset

* Option: save a histogram of the p-value distribution
"""
import h5py as h
import logging
import multiprocessing
import numpy as np

from math import ceil
from os import getpid, remove
from pyscm.binary_attributes.classifications.hdf5 import HDF5PackedAttributeClassifications
from sklearn.feature_selection.univariate_selection import chi2
from sys import argv


def compute_chi2_parallel(dataset_file, data_paths, n_cpus=-1):
    def chi2_process(proc_id, attribute_idx, tmp_file_prefix):
        logging.debug("Process " + str(proc_id) + " started.")

        # Load the attribute values for each example and the labels
        proc_dataset = h.File(dataset_file, "r")
        attribute_classifications = HDF5PackedAttributeClassifications(
            datasets=[proc_dataset[path + "/attribute_classifications"] for path in data_paths],
            n_rows=[proc_dataset[path + "/labels"].shape[0] for path in data_paths])
        labels = [proc_dataset[path + "/labels"][...] for path in data_paths]
        labels = np.asarray([item for sublist in labels for item in sublist], dtype=np.uint8)

        # Create the output HDF5 dataset
        proc_output_file = h.File(tmp_file_prefix + "%d_to_%d.tmp" % (min(attribute_idx), max(attribute_idx)), "w")
        p_values = proc_output_file.create_dataset("p_values", shape=(len(attribute_idx),), dtype=np.float64)

        # Compute the p-values in blocks to limit memory usage
        memory_block_size = 10000
        n_memory_blocks = int(ceil(1.0 * len(attribute_idx) / memory_block_size))
        for mem_block_idx in xrange(n_memory_blocks):
            memory_block_ac = attribute_classifications.get_columns(attribute_idx[mem_block_idx * memory_block_size : (mem_block_idx + 1) * memory_block_size])
            logging.debug("Process " + str(proc_id) + " done loading block " + str(mem_block_idx + 1) + " of " + str(
                n_memory_blocks))
            p_values[mem_block_idx * memory_block_size: (mem_block_idx + 1) * memory_block_size] = \
                chi2(memory_block_ac, labels)[1]

        proc_output_file.close()
        proc_dataset.close()


    if n_cpus == -1:
        n_cpus = multiprocessing.cpu_count()

    dataset = h.File(dataset_file, "r")
    presence_attribute_idx = np.where(dataset["type_by_attribute"][...] == 0)[0]
    n_presence_attributes = len(presence_attribute_idx)
    n_attributes = dataset["type_by_attribute"].shape[0]
    dataset.close()
    del dataset

    tmp_file_prefix = "%d_"%getpid()

    block_size = int(ceil(1.0 * n_presence_attributes / n_cpus))

    procs = []
    for i in xrange(n_cpus):
        p = multiprocessing.Process(
            target=chi2_process,
            args=(i, presence_attribute_idx[i * block_size : min((i + 1) * block_size, n_presence_attributes)], tmp_file_prefix))
        procs.append(p)
        p.start()

    for p in procs:
        p.join()

    # Save the merged p-values to the file
    p_values = np.ones(n_attributes) * 9999.0
    for block_start, block_end in ((i * block_size, min((i + 1) * block_size, n_presence_attributes)) for i in xrange(n_cpus)):
        tmp_file_name = tmp_file_prefix + "%d_to_%d.tmp" % (presence_attribute_idx[block_start],
                                                            presence_attribute_idx[block_end - 1])
        tmp_file = h.File(tmp_file_name, "r")
        p_values[presence_attribute_idx[block_start:block_end]] = tmp_file["p_values"][...]
        tmp_file.close()
        remove(tmp_file_name)

    # Replace the nan values with 1.0 (I think that this happens when a kmer was in all the examples)
    p_values[np.isnan(p_values)] = 1.0

    return p_values



if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s.%(msecs)d %(levelname)s - [%(process)d] - %(module)s - %(funcName)s: %(message)s")

    dataset_file = argv[1]
    n_cpus = int(argv[2])

    logging.debug("Computing the p-values for the dataset %s" % dataset_file)

    dataset = h.File(dataset_file, "r")
    folds = [k for k in dataset["train/cv"].iterkeys() if "fold" in k]
    dataset.close()

    # Compute the p-values for the cross-validation folds
    logging.debug("Computing the p-values for the cross-validation folds")
    for fold in folds:
        logging.debug("Fold: %s" % fold)
        other_folds = [other for other in folds if other != fold]
        data_paths = ["train/cv/" + other for other in other_folds]
        p_values = compute_chi2_parallel(dataset_file, data_paths, n_cpus)
        dataset = h.File(dataset_file)
        dataset["train/cv/" + fold].create_dataset("p_values", data=p_values)
        dataset.close()

    logging.debug("Computing the p-values for the training set")
    p_values = compute_chi2_parallel(dataset_file, ["train"], n_cpus)
    dataset = h.File(dataset_file)
    dataset["train"].create_dataset("p_values", data=p_values)
    dataset.close()

    logging.debug("Completed.")