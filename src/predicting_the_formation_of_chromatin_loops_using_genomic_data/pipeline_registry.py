"""Project pipelines."""
from typing import Dict

from predicting_the_formation_of_chromatin_loops_using_genomic_data.pipelines import data_downloading
from predicting_the_formation_of_chromatin_loops_using_genomic_data.pipelines import data_preprocessing
from kedro.pipeline import Pipeline


def register_pipelines() -> Dict[str, Pipeline]:
    """Register the project's pipelines.

    Returns:
        A mapping from pipeline names to ``Pipeline`` objects.
    """
    data_downloading_pipeline = data_downloading.create_pipeline()
    data_preprocessing_pipeline = data_preprocessing.create_pipeline()

    return {
        "__default__": data_downloading_pipeline + data_preprocessing_pipeline,
        "data_downloading": data_downloading_pipeline,
        "data_preprocessing": data_preprocessing_pipeline
    }