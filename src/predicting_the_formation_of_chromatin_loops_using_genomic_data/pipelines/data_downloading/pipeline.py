"""
This is a boilerplate pipeline 'data_downloading'
generated using Kedro 0.18.6
"""

from kedro.pipeline import Pipeline, node, pipeline
from .nodes import download_data

def create_pipeline(**kwargs) -> Pipeline:
    return pipeline(
        [
            node(
                func=download_data,
                inputs="params:to_download",
                outputs=None,
                name="download_data_node",
            ),
        ]
    )
