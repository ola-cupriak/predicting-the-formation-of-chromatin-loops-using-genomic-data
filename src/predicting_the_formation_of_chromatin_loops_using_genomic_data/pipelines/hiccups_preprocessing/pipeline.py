"""
This is a boilerplate pipeline 'hiccups_preprocessing'
generated using Kedro 0.18.6
"""

from kedro.pipeline import Pipeline, node, pipeline
from .nodes import preprocess_hiccups, add_labels

def create_pipeline(**kwargs) -> Pipeline:
    return pipeline(
        [
            node(
                func=preprocess_hiccups,
                inputs="HiCCUPS_loops_anotations",
                outputs="preprocessed_HiCCUPS_loops_anotations",
                name="preprocess_HiCCUPS_loops_anotations_node",
            ),
            node(
                func=add_labels,
                inputs="preprocessed_HiCCUPS_loops_anotations",
                outputs="concat_label_HiCCUPS_loops_anotations",
                name="concat_label_HiCCUPS_loops_anotations_node",
            ),
        ]
    )
