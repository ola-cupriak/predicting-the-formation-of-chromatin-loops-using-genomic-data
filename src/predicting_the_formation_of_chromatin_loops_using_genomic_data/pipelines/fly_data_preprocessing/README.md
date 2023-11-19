# Pipeline fly_data_preprocessing

## Overview

A pipeline used to prepare fly fruit datasets for training machine learning models based on data downloaded via the data_downloading pipeline. Pipeline can be run in 3 versions, depending on which method of generating negative examples is to be used:
- data_preprocessing_n2 - random pairing of regions constituting positive examples in the same cell type,
- data_preprocessing_n3 - random pairing of open chromatin regions in the same cell type,
- data_preprocessing_n4 - positive examples are long-range chromatin loops, and negative examples are short-range chromatin loops.

# Configuration

The pipeline parameters can be passed through two configuration files: conf/base/parameters.yml and conf/base/parameters/fly_data_preprocessing.yml

The below lists describe the parameters that can be changed to obtain the desired results. Parameters that are not in the lists below should not be changed, as they ensure the correct operation of the pipeline. 

1. Parameters available in conf/base/parameters.yml:
    - fly_neg_pos_ratio - the ratio of the number of negative examples to the number of positive examples to be used when generating negative examples for fly fruit datasets

2. Parameters available in conf/base/parameters/data_preprocessing.yml:
    - fly_radius - value to be subtracted and added to the center of the anchor, in order to generate anchors of the same size (the size of the created anchors will be 2*radius),
    - fly_resolution - the resolution to be used to generate statistics for continuous functional genomics data (bigWig),
    - fly_short_long_limit - length of chromatin loop above which loops are considered to be long-range loops (ie. positive examples in the method of generating negatives n4)

## Output

The final file, containing the prepared data for all the indicated human cell types and selected negative sampling type (n1/n2/n3): 04_feature/D_melanogaster/negative_sampling_type/negative_sampling_type_concatenated_combined_functional_genomics_data.parquet

## Example

Running only the data human data preprocessing pipeline with negative sampling method n1:

    kedro run -p fly_data_preprocessing_n2

Final file will be saved as: 04_feature/D_melanogaster/n2/n2_concatenated_combined_functional_genomics_data.parquet


