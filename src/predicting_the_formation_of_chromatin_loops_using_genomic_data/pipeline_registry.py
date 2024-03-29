"""Project pipelines."""
from typing import Dict

from predicting_the_formation_of_chromatin_loops_using_genomic_data.pipelines import (
    data_downloading,
)
from predicting_the_formation_of_chromatin_loops_using_genomic_data.pipelines import (
    data_preprocessing,
)
from predicting_the_formation_of_chromatin_loops_using_genomic_data.pipelines import (
    fly_data_preprocessing,
)
from predicting_the_formation_of_chromatin_loops_using_genomic_data.pipelines import (
    model_training,
)
from predicting_the_formation_of_chromatin_loops_using_genomic_data.pipelines import (
    fly_model_training,
)

from kedro.pipeline import Pipeline


def register_pipelines() -> Dict[str, Pipeline]:
    """Register the project's pipelines.

    Returns:
        A mapping from pipeline names to ``Pipeline`` objects.
    """
    data_downloading_pipeline = data_downloading.create_pipeline()
    data_preprocessing_n1_pipeline = data_preprocessing.create_pipeline(
        neg_sampling_type="n1"
    )
    data_preprocessing_n2_pipeline = data_preprocessing.create_pipeline(
        neg_sampling_type="n2"
    )
    data_preprocessing_n3_pipeline = data_preprocessing.create_pipeline(
        neg_sampling_type="n3"
    )

    fly_data_preprocessing_n2_pipeline = fly_data_preprocessing.create_pipeline(
        neg_sampling_type="n2"
    )
    fly_data_preprocessing_n3_pipeline = fly_data_preprocessing.create_pipeline(
        neg_sampling_type="n3"
    )
    fly_data_preprocessing_n4_pipeline = fly_data_preprocessing.create_pipeline(
        neg_sampling_type="n4"
    )

    model_training_within_pipeline = model_training.create_pipeline(mtype="within")
    model_training_across_pipeline = model_training.create_pipeline(mtype="across")

    fly_model_training_pipeline = fly_model_training.create_pipeline()

    return {
        "__default__": data_downloading_pipeline
        + data_preprocessing_n1_pipeline
        + model_training_within_pipeline,
        "data_downloading": data_downloading_pipeline,
        "data_preprocessing_n1": data_preprocessing_n1_pipeline,
        "data_preprocessing_n2": data_preprocessing_n2_pipeline,
        "data_preprocessing_n3": data_preprocessing_n3_pipeline,
        "fly_data_preprocessing_n2": fly_data_preprocessing_n2_pipeline,
        "fly_data_preprocessing_n3": fly_data_preprocessing_n3_pipeline,
        "fly_data_preprocessing_n4": fly_data_preprocessing_n4_pipeline,
        "model_training_across_cells": model_training_across_pipeline,
        "model_training_within_cells": model_training_within_pipeline,
        "fly_model_training": fly_model_training_pipeline,
    }
