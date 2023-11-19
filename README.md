# Predicting the formation of chromatin loops using genomic data

## Contents
  * [General Description](#general-description)
  * [Installation](#installation)
  * [General Information](#general-information)
  * [Configuration](#configuration)
  * [Examples](#examples)

## General Description

The project aims to use machine learning models to predict chromatin loop formation based on genomic data. The models  are designed to perform binary classification to assess whether two genomic regions form a chromatin loop or not, based on DNase-seq, CTCF ChIP-seq and vertebrate motifs data from 7 human cell types or scATAC-seq, various ChIP-seq and insect motifs data from from the central nervous system of the fruit fly. The project consists of 5 pipelines: for data downloading, for human data preprocessing, for fruit fly data preprocessing, for training and evaluation of models for human data and for training and evaluation of models on fruit fly data.

## Requirments
  * [Miniconda](https://docs.conda.io/en/latest/miniconda.html)

## Installation
The best way to install the project with all necessary dependencies is to use conda and the environment.yml file.

    conda env create -f environment.yml
    
## General Information
1. The project was built using Kedro (https://kedro.org/) - Python framework for creating reproducible, maintainable and modular data-science code.

2. Detailed documentation and the main code of each pipeline can be found in the node.py and pipeline.py files located in the directory of each of the pipelines found in the directory "src/predicting_the_formation_of_chromatin_loops_using_genomic_data/pipelines/".

3. There are 3 methods of generating negative examples during data preprocessing:
 - anchors_from_other_cell_types (n1) - random pairing of regions constituting positive examples in other cell types,
 - new_pairs_of_anchors (n2) - random pairing of regions constituting positive examples in the same cell type,
 - open_chromatin_regions (n3) - random pairing of open chromatin regions in the same cell type,
 - FLY_short_loops_as_negatives (n4) - positive examples are long-range chromatin loops, and negative examples are short-range chromatin loops. 

4. There aretypes of machine learning models:
    - within - models trained and evaluated on data from the same cell type,
    - across - models trained on data from multiple cell types and evaluated on data from a new cell type unknown to the model.

5. Available algorithms for machine learning models:
    - logistic regression,
    - decision tree,
    - random forest,
    - LightGBM.

6. 5 pipelines were created for various tasks. Some of them can be run in different versions (considering point 3. and 4.):
    - data_downloading - pipeline for downloading data specified in conf/base/parameters/data_downloading.yml file,
    - data_preprocessing - pipeline for preparing human data. It can be run in 3 versions, depending on which method of generating negative examples is to be used:
        - data_preprocessing_n1,
        - data_preprocessing_n2,
        - data_preprocessing_n3.
    - fly_data_preprocessing - pipeline for preparing fruit fly data. It can be run in 3 versions, depending on which method of generating negative examples is to be used.
        - fly_data_preprocessing_n2,
        - fly_data_preprocessing_n3,
        - fly_data_preprocessing_n4.
    - model_training - pipeline for training and evaluating machine learning models on human data. It can be run in 2 versions, depending on whether to train within or across models:
        - model_training_within_cells,
        - model_training_across_cells.
    - fly_model_training -  pipeline for training and evaluating machine learning models (type: within) on fruit fly data.

7. Various configurations of pipelining launches are possible. Do not combine pipelines for fruit fly (with the prefix "fly") with pipelines for humans (without the prefix "fly"). Only the data_downloading pipeline is universal for both organisms. By default, a configuration consisting of the pipelines: data_downloading + data_preprocessing_n1 + model_training_within_cells.

8. Tracking of models is possible with MLflow (local host).

9. A summary of the metrics (accuracy, auc) and learning times is saved to the file run_stats.txt in the data/08_reportiong directory.

10. The exact values of the metrics, optimalization info and feature importances for each trained model are saved in the directory corresponding to the experiment in the mlruns directory.

## Configuration

Parameter configurations for each pipeline are discussed in the corresponding README files.

## Examples

1. Run with the default settings (data_downloading + data_preprocessing_n1 + model_training_within_cells):

 kedro run
 
3. Run only 1 pipeline (e.g., preparation of human data with n1 negative generation method): kedro run -p data_preprocessing_n1
4. Enable MLflow UI: kedro mlflow ui

