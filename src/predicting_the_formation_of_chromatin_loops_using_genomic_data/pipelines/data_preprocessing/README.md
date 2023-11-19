# Pipeline data_preprocessing

## Overview

A pipeline used to prepare human datasets for training machine learning models based on data downloaded via the data_downloading pipeline. Pipeline can be run in 3 versions, depending on which method of generating negative examples is to be used:

    - data_preprocessing_n1 - random pairing of regions constituting positive examples in other cell types,
    - data_preprocessing_n2 - random pairing of regions constituting positive examples in the same cell type,
    - data_preprocessing_n3 - random pairing of open chromatin regions in the same cell type.

# Configuration

The pipeline parameters can be passed through two configuration files: conf/base/parameters.yml and conf/base/parameters/data_preprocessing.yml

The below lists describe the parameters that can be changed to obtain the desired results. Parameters that are not in the lists below should not be changed, as they ensure the correct operation of the pipeline. 

1. Parameters available in conf/base/parameters.yml:
    - human_neg_pos_ratio - the ratio of the number of negative examples to the number of positive examples to be used when generating negative examples for human datasets

2. Parameters available in conf/base/parameters/data_preprocessing.yml:
    - radius - value to be subtracted and added to the center of the anchor, in order to generate anchors of the same size (the size of the created anchors will be 2*radius),
    - resolution - the resolution to be used to generate statistics for continuous functional genomics data (bigWig),
    - cell_types_to_use - list of cell types to be included in the data.

## Output

The final file, containing the prepared data for all the indicated human cell types and selected negative sampling type (n1/n2/n3): 04_feature/Homo_sapiens/negative_sampling_type/negative_sampling_type_concatenated_combined_functional_genomics_data.parquet

## Example

Running only the data human data preprocessing pipeline with negative sampling method n1:

    kedro run -p data_preprocessing_n1

Final file will be saved as: 04_feature/Homo_sapiens/n1/n1_concatenated_combined_functional_genomics_data.parquet


