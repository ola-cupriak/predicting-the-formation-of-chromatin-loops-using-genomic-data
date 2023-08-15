"""
This is a boilerplate pipeline 'data_downloading'
generated using Kedro 0.18.6
"""

from kedro.pipeline import Pipeline, node, pipeline
from .nodes import (
    download_data,
    simplify_human_genome_file,
    simplify_flies_genome_file,
    map_dm3_to_dm6,
)


def create_pipeline(**kwargs) -> Pipeline:
    return pipeline(
        [
            node(
                func=download_data,
                inputs="params:to_download",
                outputs=[
                    "path_hg19_raw",
                    "path_dm6_raw",
                    "cells2names_param",
                    "FLY_cells2names_param",
                    "dm3_to_dm6_chain",
                ],
                name="download_data_node",
            ),
            node(
                func=simplify_human_genome_file,
                inputs="path_hg19_raw",
                outputs="hg19_simplified",
                name="simplify_hg19_node",
            ),
            node(
                func=simplify_flies_genome_file,
                inputs="path_dm6_raw",
                outputs="dm6_simplified",
                name="simplify_dm6_node",
            ),
            node(
                func=map_dm3_to_dm6,
                inputs=[
                    "FLY_cells2names_param",
                    "params:dm3todm6_mapping",
                    "dm3_to_dm6_chain",
                ],
                outputs=None,
                name="map_dm3_to_dm6_node",
            ),
        ]
    )
