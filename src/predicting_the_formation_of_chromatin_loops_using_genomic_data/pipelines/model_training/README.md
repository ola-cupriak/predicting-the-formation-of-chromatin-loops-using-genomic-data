# Pipeline model_training

## Overview

Pipeline used to train and evaluate machine learning models on human datasets. Pipeline can be run in 2 versions:
- model_training_within_cells - models trained and evaluated on data from the same cell type,
- model_training_across_cells - models trained on data from multiple cell types and evaluated on data from a new cell type unknown to the model.

## Configuration
The pipeline parameters can be passed through configuration file: conf/base/parameters/model_training.yml

The below lists describe the parameters that can be changed to obtain the desired results. Parameters that are not in the lists below should not be changed, as they ensure the correct operation of the pipeline.

- run_name: name of experiment to be recorded in MLflow
- random_state
- data_fraction: fraction of data to be used for training and evaluation of models (the same percentage of positive and negative examples will be drawn); float in the range 0-1
- cell_types: list of cell types to be used for training and evaluation of models
- features_include_only: the list of features (columns in the data) to be used for training and evaluation of models.
- features_exclude: the list of features (columns in the data) to be excluded from training and evaluation of models.

IMPORTANT: If both lists, features_include_only and features_exclude, are empty, all features available in the data will be used for training and evaluation of models.

Parameters for individual models:
- run: decides whether the model will be launched; bool
- optimize: decides whether the optimization of the model's hyperparameters will be launched; bool
- optim_time: maximum optimization time, in seconds; int
- n_trials: maximum number of optimization trials; int
- params: dictionary of individual parameters to optimize for a given model

You can read about the parameters for models in scikit-learn and Pytorch-tabular documentation.

## Pipeline inputs

<!---
The list of pipeline inputs.
-->

## Pipeline outputs

<!---
The list of pipeline outputs.
-->

