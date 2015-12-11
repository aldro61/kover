---
title: Kover learning engine
tags:
keywords: learning, engine, set covering machine, SCM, kover, genomics, k-mer, machine learning
last_updated: December 11, 2015
summary: "Overview of the Kover learning engine"
---

## Learning models

```
usage: kover learn [-h] --dataset DATASET --split SPLIT --model-type
                   {conjunction,disjunction} [{conjunction,disjunction} ...]
                   --p P [P ...] --max-rules MAX_RULES
                   [--hp-choice {bound,cv,none}] [--n-cpu N_CPU]
                   [--output-dir OUTPUT_DIR] [-x] [-v]

Learn a model from data

optional arguments:
  -h, --help            show this help message and exit
  --dataset DATASET     The Kover dataset to learn from.
  --split SPLIT         The identifier of the split of the dataset that must
                        be learnt from.
  --model-type {conjunction,disjunction} [{conjunction,disjunction} ...]
                        Hyperparameter: The type of model (conjunction or
                        disjunction) to learn. Single value or multiple space
                        separated values.
  --p P [P ...]         Hyperparameter: The value of the trade-off used to
                        score the rules. Single value or multiple space
                        separated values.
  --max-rules MAX_RULES
                        The maximum number of rules to include in a model.
  --hp-choice {bound,cv,none}
                        The strategy used to select the hyperparameter values.
  --n-cpu N_CPU         The number of CPUs used to select the hyperparameter
                        values. Make sure your computer has enough RAM and
                        that your storage device is not a bottleneck (see
                        documentation).
  --output-dir OUTPUT_DIR
                        The directory in which to store Kover's output. Will
                        be created if it does not exist.
  -x, --progress        Shows a progress bar for the execution.
  -v, --verbose         Sets the verbosity level.
```