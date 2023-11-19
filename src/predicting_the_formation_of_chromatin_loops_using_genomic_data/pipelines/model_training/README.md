# Pipeline model_training

## Overview

Pipeline used to train and evaluate machine learning models on human datasets. Pipeline can be run in 2 versions:
- model_training_within_cells - models trained and evaluated on data from the same cell type,
- model_training_across_cells - models trained on data from multiple cell types and evaluated on data from a new cell type unknown to the model.

The best way to track the results received is to use the MLFlow UI. You can run it with the command: 
``` kedro mlflow ui ```

## Configuration
The pipeline parameters can be passed through configuration file: conf/base/parameters/model_training.yml

The below lists describe the parameters that can be changed to obtain the desired results. Parameters that are not in the lists below should not be changed, as they ensure the correct operation of the pipeline.

1. General parameters:
- run_name - name of experiment to be recorded in MLflow; it is **good practice** to set unique names for experiments that **start with "within" or "across" word** depending on the type of models being trained,
- data_fraction - fraction of data to be used for training and evaluation of models (the same percentage of positive and negative examples will be drawn); float in the range 0-1,
- neg_sampling_type - type of data to be used (depending on the method of generating negatives); choice: n1/n2/n3,
- test_size - size of the test dataset; float in the range 0-1,
- validation_size - size of the validation dataset; used only when hyperparameter optimization is running; float in the range 0-1,
- cell_types - list of cell types to be used for training and evaluation of models,
- features_include_only: the list of features (columns in the data) to be used for training and evaluation of models,
- features_exclude: the list of features (columns in the data) to be excluded from training and evaluation of models.

If both lists, features_include_only and features_exclude, are empty, all features available in the data will be used for training and evaluation of models.
If both lists, features_include_only and features_exclude, are not empty, only features_include_only will be taken into account.

2. Parameters for individual models:
- run - decides whether the model will be launched; bool
- optimize - decides whether the optimization of the model's hyperparameters will be launched; bool
- optim_time - maximum optimization time, in seconds; int
- n_trials - maximum number of optimization trials; int
- params:
    - if optimize=True - dictionary of individual parameters to optimize for a given model (each parameter passed should have a specified type - int/float/categorical),
    - if optimize=False - should be left empty or contain specific parameter values to be used (e.g. max_features: 10)

For a detailed description of the parameters that can be passed, see the scikit-learn documentation for specific model.

