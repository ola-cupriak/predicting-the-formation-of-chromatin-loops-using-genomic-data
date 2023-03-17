"""
This is a boilerplate pipeline 'hiccups_preprocessing'
generated using Kedro 0.18.6
"""

from kedro.pipeline import Pipeline, node, pipeline
from .nodes import read_hiccups, add_labels, read_signals, count_signals

def create_pipeline(**kwargs) -> Pipeline:
    return pipeline(
        [
            node(
                func=read_hiccups,
                inputs="HiCCUPS_loops_anotations",
                outputs="readed_HiCCUPS_loops_anotations",
                name="read_HiCCUPS_loops_anotations_node",
            ),
            node(
                func=add_labels,
                inputs="readed_HiCCUPS_loops_anotations",
                outputs="concat_label_HiCCUPS_loops_anotations",
                name="concat_label_HiCCUPS_loops_anotations_node",
            ),
            node(
                func=read_signals,
                inputs="DNAse_seq_signals",
                outputs="readed_DNase_seq_signals",
                name="read_DNase_seq_signals_node",
            ),
            # node(
            #     func=read_signals,
            #     inputs="CTCF_ChIP_seq_signals",
            #     outputs="readed_CTCF_ChIP_seq_signals",
            #     name="read_CTCF_ChIP_seq_signals_node",
            # ),
            node(
                func=count_signals,
                inputs=["concat_label_HiCCUPS_loops_anotations", "readed_DNase_seq_signals", "params:DNAse-seq_name", "params:radius"],
                outputs="combined_functional_genomics_data",
                name="combine_functional_genomics_data_node",
            ),
        ]
    )
