"""
This is a boilerplate pipeline 'data_downloading'
generated using Kedro 0.18.6
"""

from kedro.pipeline import Pipeline, node, pipeline
from .nodes import download_data, simplify_genome_file

def create_pipeline(**kwargs) -> Pipeline:
    return pipeline(
        [
            node(
                func=download_data,
                inputs="params:to_download",
                outputs="path_hg19_raw",
                name="download_data_node",
            ),
            node(
                func=simplify_genome_file,
                inputs="path_hg19_raw",
                outputs="hg19_simplified",
                name="simplify_hg19_node",
            ),
        ]
    )
