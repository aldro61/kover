#!/usr/bin/env python
import logging
from math import ceil, floor

import h5py as h
import numpy as np
from numpy import uint8


def get_row_2(packed_data, inner_row_idx):
    if packed_data.dtype == np.uint8:
        pack_size = 8
        pack_type = np.uint8
    elif packed_data.dtype == np.uint16:
        pack_size = 16
        pack_type = np.uint16
    elif packed_data.dtype == np.uint32:
        pack_size = 32
        pack_type = np.uint32
    elif packed_data.dtype == np.uint64:
        pack_size = 64
        pack_type = np.uint64
    else:
        raise ValueError("Supported data types are 8-bit, 16-bit, 32-bit and 64-bit integers.")

    return np.asarray(np.bitwise_and(np.right_shift(packed_data, pack_type(pack_size - inner_row_idx - 1)), 1),
                      dtype=np.uint8)

def set_row(packed_data, unpacked_row, unpacked_row_idx, clear_before_set=True):
    if packed_data.dtype == np.uint8:
        pack_size = 8
        pack_type = np.uint8
    elif packed_data.dtype == np.uint16:
        pack_size = 16
        pack_type = np.uint16
    elif packed_data.dtype == np.uint32:
        pack_size = 32
        pack_type = np.uint32
    elif packed_data.dtype == np.uint64:
        pack_size = 64
        pack_type = np.uint64
    else:
        raise ValueError("Supported data types are 8-bit, 16-bit, 32-bit and 64-bit integers.")

    packed_data_row = int(floor(1.0 * unpacked_row_idx / pack_size))

    max_idx = packed_data.shape[0] * pack_size - 1
    if unpacked_row_idx > max_idx:
        raise IndexError("Index %d is out of bounds for bit array of size %d." % (unpacked_row_idx, max_idx + 1))

    tmp = np.asarray(unpacked_row, dtype=pack_type)
    tmp = np.left_shift(tmp, pack_type(pack_size - (unpacked_row_idx - (packed_data_row * pack_size)) - 1))

    if clear_before_set:
        mask = np.ones(packed_data.shape[1], dtype=pack_type)
        mask = np.left_shift(mask, pack_type(pack_size - (unpacked_row_idx - (packed_data_row * pack_size)) - 1))
        np.bitwise_not(mask, out=mask)
        np.bitwise_and(packed_data[packed_data_row], mask, out=packed_data[packed_data_row])

    np.bitwise_or(packed_data[packed_data_row], tmp, out=packed_data[packed_data_row])
    return packed_data


def split_into_folds(fold_prefix, n_folds, destination, attribute_classifications, labels, example_identifiers,
                     random_generator, gzip, chunk_size):
    """
    We store only the test data for each fold, since it can be combined to obtain the fold training data.
    """
    # Note: the examples are already sorted in label order to optimize training speed.

    if attribute_classifications.dtype == np.uint32:
        pack_size = 32
        pack_dtype = np.uint32
    elif attribute_classifications.dtype == np.uint64:
        pack_size = 64
        pack_dtype = np.uint64
    else:
        raise ValueError("Supported data types for packing are 8-bit, 16-bit, 32-bit and 64-bit integers.")

    n_examples = len(labels)
    n_attributes = attribute_classifications.shape[1]

    # Randomly assign a fold to each example
    fold_by_example = np.arange(n_examples) % n_folds
    random_generator.shuffle(fold_by_example)

    # Create the group for each fold
    fold_groups = [destination.create_group(fold_prefix + str(i + 1)) for i in xrange(n_folds)]

    # Fold data buffers for attribute classifications
    fold_ac_datasets = [fold_groups[i].create_dataset("attribute_classifications",
                                                      shape=(int(ceil(
                                                          1.0 * len(np.where(fold_by_example == i)[0]) / pack_size)),
                                                             n_attributes),
                                                      dtype=pack_dtype,
                                                      compression="gzip" if gzip > 0 else None,
                                                      compression_opts=gzip if gzip > 0 else None,
                                                      chunks=chunk_size) for i in xrange(n_folds)]
    fold_ac_buffers = np.zeros((n_folds, n_attributes), dtype=pack_dtype)
    fold_ac_buffer_packed_rows = [0] * n_folds
    fold_output_ac_current_row = [0] * n_folds

    # Fold label buffers (Saved unpacked at the end only)
    fold_label_buffers = [[] for _ in xrange(n_folds)]
    fold_label_datasets = [
        fold_groups[i].create_dataset("labels", shape=(len(np.where(fold_by_example == i)[0]),), dtype=np.uint8) for i
        in xrange(n_folds)]

    # Fold example identifiers (indexes with respect to the master list that is stored at the root of the file)
    fold_example_identifier_buffers = [[] for _ in xrange(n_folds)]
    fold_example_identifier_datasets = [
        fold_groups[i].create_dataset("example_identifiers", shape=(len(np.where(fold_by_example == i)[0]),),
                                      dtype=example_identifiers.dtype) for i in xrange(n_folds)]

    example_count = 0
    for packed_block in attribute_classifications:
        for packed_idx in xrange(pack_size):
            if example_count == n_examples:
                logging.debug("Stopping since all examples have been packed.")
                break  # If all examples have been packed, stop

            fold = fold_by_example[example_count]

            logging.debug("Sample %d/%d - Destination: Fold #%d - Fold buffer status: %d/%d" %
                          (example_count + 1, n_examples, fold + 1, fold_ac_buffer_packed_rows[fold], pack_size))

            unpacked_row = get_row_2(packed_block, packed_idx)

            fold_ac_buffers[fold] = set_row(fold_ac_buffers[fold].reshape(1, -1),
                                            unpacked_row,
                                            fold_ac_buffer_packed_rows[fold],
                                            False)
            fold_ac_buffer_packed_rows[fold] += 1

            fold_label_buffers[fold].append(labels[example_count])
            fold_example_identifier_buffers[fold].append(example_identifiers[example_count])

            if fold_ac_buffer_packed_rows[fold] == pack_size:
                logging.debug("Flushing buffer for fold %d" % (fold + 1))
                # Flush buffer
                fold_ac_datasets[fold][fold_output_ac_current_row[fold]] = fold_ac_buffers[fold]
                # Reset buffer
                fold_ac_buffers[fold] = np.zeros((1, n_attributes), dtype=pack_dtype)
                fold_ac_buffer_packed_rows[fold] = 0
                # Increment output index
                fold_output_ac_current_row[fold] += 1

            example_count += 1

    logging.debug("Flushing all buffers.")
    for fold in xrange(n_folds):
        # Flush the attribute classifications buffers
        if fold_ac_buffer_packed_rows[fold] != 0:
            fold_ac_datasets[fold][fold_output_ac_current_row[fold]] = fold_ac_buffers[fold]

        # Flush the label buffers
        fold_label_datasets[fold][...] = fold_label_buffers[fold]

        # Flush the example identifier buffers
        fold_example_identifier_datasets[fold][...] = fold_example_identifier_buffers[fold]


def split_train_test(input_file, output_file, train_prop, random_generator, gzip):
    # Note: the examples are already sorted in label order to optimize training speed.

    if input_file["attribute_classifications"].dtype == np.uint64:
        pack_size = 64
        pack_dtype = np.uint64
    elif input_file["attribute_classifications"].dtype == np.uint32:
        pack_size = 32
        pack_dtype = np.uint32
    else:
        raise RuntimeError(
            "Unsupported pack type for the attribute_classification matrix. Supported types are np.uint32 and np.uint64.")

    attribute_classifications = input_file["attribute_classifications"]
    labels = input_file["labels"][...]
    example_identifiers = np.arange(len(labels))  # Reference to the example identifier dataset
    n_attributes = input_file["attribute_classifications"].shape[1]
    n_samples = len(labels)
    n_train_samples = int(n_samples * train_prop)
    n_test_samples = n_samples - n_train_samples

    # Create train group
    logging.debug("Creating the training datasets")
    train_group = output_file.create_group("train")
    train_attribute_classifications = train_group.create_dataset("attribute_classifications",
                                                                 dtype=input_file["attribute_classifications"].dtype,
                                                                 shape=(int(ceil(1.0 * n_train_samples / pack_size)),
                                                                        n_attributes),
                                                                 chunks=input_file["attribute_classifications"].chunks,
                                                                 compression="gzip" if gzip > 0 else None,
                                                                 compression_opts=gzip if gzip > 0 else None)
    train_labels = train_group.create_dataset("labels",
                                              dtype=input_file["labels"].dtype,
                                              shape=(n_train_samples,))
    train_example_identifiers = train_group.create_dataset("example_identifiers",
                                                           dtype=np.uint32,
                                                           shape=(n_train_samples,))

    # Create test group
    logging.debug("Creating the testing datasets")
    test_group = output_file.create_group("test")
    test_attribute_classifications = test_group.create_dataset("attribute_classifications",
                                                               dtype=input_file["attribute_classifications"].dtype,
                                                               shape=(int(ceil(1.0 * n_test_samples / pack_size)),
                                                                      n_attributes),
                                                               chunks=input_file["attribute_classifications"].chunks,
                                                               compression="gzip" if gzip > 0 else None,
                                                               compression_opts=gzip if gzip > 0 else None)
    test_labels = test_group.create_dataset("labels",
                                            dtype=input_file["labels"].dtype,
                                            shape=(n_test_samples,))
    test_example_identifiers = test_group.create_dataset("example_identifiers",
                                                         dtype=np.uint32,
                                                         shape=(n_test_samples,))

    # Randomly assign a set (train or test) to each sample
    logging.debug("Randomly assigning a set (train/test) to the samples")
    sample_destination = np.zeros(n_samples, dtype=np.bool)
    idx = np.arange(len(labels))
    random_generator.shuffle(idx)
    sample_destination[idx[: n_train_samples]] = True  # True is for training and False for testing

    # Write the examples to the correct datasets
    # TODO: We could use column blocks to reduce memory usage
    logging.debug("Initializing training set buffers")
    train_ac_buffer = np.zeros((1, n_attributes), dtype=pack_dtype)
    train_ac_buffer_packed_rows = 0
    train_ac_output_current_row = 0
    train_label_buffer = []
    train_example_identifiers_buffer = []

    logging.debug("Initializing testing set buffers")
    test_ac_buffer = np.zeros((1, n_attributes), dtype=pack_dtype)
    test_ac_buffer_packed_rows = 0
    test_ac_output_current_row = 0
    test_label_buffer = []
    test_example_identifiers_buffer = []

    example_count = 0
    for packed_block in attribute_classifications:
        for packed_idx in xrange(pack_size):
            if example_count == n_samples:
                logging.debug("Stopping since all examples have been packed.")
                break  # If all examples have been packed, stop

            logging.debug("Sample %d/%d - Destination: %s set - Buffer status: train=(%d/%d), test=(%d/%d)" %
                          (example_count + 1, n_samples, "training" if sample_destination[example_count] else "testing",
                           train_ac_buffer_packed_rows, pack_size, test_ac_buffer_packed_rows, pack_size))

            unpacked_row = get_row_2(packed_block, packed_idx)

            if sample_destination[example_count]:  # Add the example to the training set
                train_ac_buffer = set_row(train_ac_buffer.reshape(1, -1),
                                          unpacked_row,
                                          train_ac_buffer_packed_rows,
                                          False)
                train_ac_buffer_packed_rows += 1
                train_label_buffer.append(labels[example_count])
                train_example_identifiers_buffer.append(example_identifiers[example_count])

                if train_ac_buffer_packed_rows == pack_size:
                    logging.debug("Flushing the training set ac buffer")
                    # Flush buffer
                    train_attribute_classifications[train_ac_output_current_row] = train_ac_buffer
                    # Reset buffer
                    train_ac_buffer = np.zeros((1, n_attributes), dtype=pack_dtype)
                    train_ac_buffer_packed_rows = 0
                    # Increment output index
                    train_ac_output_current_row += 1

            else:  # Add the example to the testing set
                test_ac_buffer = set_row(test_ac_buffer.reshape(1, -1),
                                         unpacked_row,
                                         test_ac_buffer_packed_rows,
                                         False)
                test_ac_buffer_packed_rows += 1
                test_label_buffer.append(labels[example_count])
                test_example_identifiers_buffer.append(example_identifiers[example_count])

                if test_ac_buffer_packed_rows == pack_size:
                    logging.debug("Flushing the testing set ac buffer")
                    # Flush buffer
                    test_attribute_classifications[test_ac_output_current_row] = train_ac_buffer
                    # Reset buffer
                    test_ac_buffer = np.zeros((1, n_attributes), dtype=pack_dtype)
                    test_ac_buffer_packed_rows = 0
                    # Increment output index
                    test_ac_output_current_row += 1

            example_count += 1

    # Flush all buffers if needed
    logging.debug("Flushing all buffers.")
    if train_ac_buffer_packed_rows > 0:
        train_attribute_classifications[train_ac_output_current_row] = train_ac_buffer
        train_labels[...] = train_label_buffer
        train_example_identifiers[...] = train_example_identifiers_buffer
    if test_ac_buffer_packed_rows > 0:
        test_attribute_classifications[test_ac_output_current_row] = test_ac_buffer
        test_labels[...] = test_label_buffer
        test_example_identifiers[...] = test_example_identifiers_buffer


def split(input, output, train_prop=0.5, random_seed=42, n_folds=5, gzip=4):
    logging.debug("Initializing a random number generator with seed %d" % random_seed)
    random_generator = np.random.RandomState(seed=random_seed)

    # Create the output file
    logging.debug("Creating the output file.")
    output_file = h.File(output, "w")

    # Copy the general stuff
    input_file = h.File(input, "r")

    logging.debug("Copying the kmer list to %s", output_file.filename)
    output_file.copy(input_file["kmers"], "kmers")

    logging.debug("Copying the example identifiers to %s", output_file.filename)
    output_file.copy(input_file["example_identifiers"], "example_identifiers")

    logging.debug("Copying attribute types to %s", output_file.filename)
    output_file.copy(input_file["type_by_attribute"], "type_by_attribute")

    logging.debug("Copying attribute kmers to %s", output_file.filename)
    output_file.copy(input_file["kmer_by_attribute"], "kmer_by_attribute")

    # Split the training and testing sets
    logging.debug("Splitting the training (prop=%.2f) and testing (prop=%.2f) sets with random seed %d." % (
    train_prop, 1.0 - train_prop, random_seed))
    split_train_test(input_file=input_file,
                     output_file=output_file,
                     train_prop=train_prop,
                     random_generator=random_generator,
                     gzip=gzip)

    if n_folds > 1:
        # Split the training set into random folds using the training attribute classifications
        logging.debug(
            "Splitting the training set into %d random non-overlapping cross-validation folds with random seed %d." % (
            n_folds, random_seed))
        train_group = output_file["train"]
        cv_group = train_group.create_group("cv")
        split_into_folds(fold_prefix="fold_",
                         n_folds=n_folds,
                         destination=cv_group,
                         attribute_classifications=train_group["attribute_classifications"],
                         labels=train_group["labels"][...],
                         example_identifiers=train_group["example_identifiers"],
                         chunk_size=train_group["attribute_classifications"].chunks,
                         gzip=gzip,
                         random_generator=random_generator)

    logging.debug("Completed.")