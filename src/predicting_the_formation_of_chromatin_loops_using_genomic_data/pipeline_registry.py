"""Project pipelines."""
from typing import Dict

from predicting_the_formation_of_chromatin_loops_using_genomic_data.pipelines import data_downloading
from predicting_the_formation_of_chromatin_loops_using_genomic_data.pipelines import data_preprocessing
from predicting_the_formation_of_chromatin_loops_using_genomic_data.pipelines import model_training
from kedro.pipeline import Pipeline


def register_pipelines() -> Dict[str, Pipeline]:
    """Register the project's pipelines.

    Returns:
        A mapping from pipeline names to ``Pipeline`` objects.
    """
    data_downloading_pipeline = data_downloading.create_pipeline()
    data_preprocessing_pipeline = data_preprocessing.create_pipeline()
    # model_training_GM12878_pipeline = model_training.create_pipeline(cell_type='GM12878')
    # model_training_K562_pipeline = model_training.create_pipeline(cell_type='K562')
    # model_training_HMEC_pipeline = model_training.create_pipeline(cell_type='HMEC')
    # model_training_HUVEC_pipeline = model_training.create_pipeline(cell_type='HUVEC')
    # model_training_NHEK_pipeline = model_training.create_pipeline(cell_type='NHEK')
    # model_training_HeLa_pipeline = model_training.create_pipeline(cell_type='HeLa')
    # model_training_IMR90_pipeline = model_training.create_pipeline(cell_type='IMR90')
    model_training_within_pipeline = model_training.create_pipeline(type='within')
    model_training_across_pipeline = model_training.create_pipeline(type='across')

    return {
        "__default__": data_downloading_pipeline + data_preprocessing_pipeline + model_training_across_pipeline,
        "data_downloading": data_downloading_pipeline,
        "data_preprocessing": data_preprocessing_pipeline,
        # "model_training_within_cell": (
        #     model_training_GM12878_pipeline
        #     + model_training_K562_pipeline
        #     + model_training_HMEC_pipeline
        #     + model_training_HUVEC_pipeline
        #     + model_training_NHEK_pipeline
        #     + model_training_HeLa_pipeline
        #     + model_training_IMR90_pipeline
        # ),
        "model_training_across_cells": model_training_across_pipeline,
        "model_training_within_cell": model_training_within_pipeline,
    }