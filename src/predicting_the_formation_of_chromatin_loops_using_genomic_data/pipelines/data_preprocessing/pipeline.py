"""
This is a boilerplate pipeline 'data_preprocessing'
generated using Kedro 0.18.6
"""

from kedro.pipeline import Pipeline, node, pipeline
from .nodes import read_hic, add_labels, read_peaks, count_peaks, read_bigWig, add_bigWig_data, find_all_motifs

def create_pipeline(**kwargs) -> Pipeline:
    return pipeline(
        [
            node(
                func=read_hic,
                inputs=["HiC_loops_annoatations", "cells2names", "params:HiC_data"],
                outputs="readed_HiC_loops_anotations",
                name="read_HiC_loops_anotations_node",
            ),
            node(
                func=add_labels,
                inputs="readed_HiC_loops_anotations",
                outputs="concat_label_HiC_loops_anotations",
                name="concat_label_HiC_loops_anotations_node",
            ),
            node(
                func=read_peaks,
                inputs=["DNAse_seq_peaks", "cells2names", "params:DNase-seq_peaks"],
                outputs="readed_DNase_seq_peaks",
                name="read_DNase_seq_peaks_node",
            ),
            # node(
            #     func=read_peaks,
            #     inputs=["CTCF_ChIP_seq_peaks", "cells2names", "params:CTCF_ChIP-seq_peaks"],
            #     outputs="readed_CTCF_ChIP_seq_peaks",
            #     name="read_CTCF_ChIP_seq_peaks_node",
            # ),
            # node(
            #     func=read_bigWig,
            #     inputs=["DNAse_seq_bigWig", "cells2names", "params:DNase-seq_bigWig"],
            #     outputs="readed_DNase_seq_bigWig",
            #     name="read_DNase_seq_bigWig_data_node",
            # ),node(
            #     func=read_bigWig,
            #     inputs=["CTCF_ChIP_seq_bigWig", "cells2names", "params:CTCF_ChIP-seq_bigWig"],
            #     outputs="readed_CTCF_ChIP_seq_bigWig",
            #     name="read_CTCF_ChIP_seq_bigWig_data_node",
            # ),
            # node(
            #     func=count_peaks,
            #     inputs=["concat_label_HiC_loops_anotations", "readed_CTCF_ChIP_seq_peaks", "params:CTCF_ChIP-seq_peaks", "params:radius"],
            #     outputs="HiC_loops_anotations_with_CTCF_ChIP_seq_peaks",
            #     name="add_CTCF_ChIP_seq_peaks_node",
            # ),
            # node(
            #     func=count_peaks,
            #     inputs=["HiC_loops_anotations_with_CTCF_ChIP_seq_peaks", "readed_DNase_seq_peaks", "params:DNase-seq_peaks", "params:radius"],
            #     outputs="HiC_loops_anotations_with_DNase_seq_peaks",
            #     name="add_DNase_seq_peaks_node",
            # ),
            # node(
            #     func=add_bigWig_data,
            #     inputs=["HiC_loops_anotations_with_DNase_seq_peaks", "readed_DNase_seq_bigWig", "params:DNase-seq_bigWig", "params:radius"],
            #     outputs="HiC_loops_anotations_with_DNase_seq_bigWig_data",
            #     name="add_DNase_seq_bigWig_data_node",
            # ),
            # node(
            #     func=add_bigWig_data,
            #     inputs=["HiC_loops_anotations_with_DNase_seq_bigWig_data", "readed_CTCF_ChIP_seq_bigWig", "params:CTCF_ChIP-seq_bigWig", "params:radius"],
            #     outputs="combined_functional_genomics_data",
            #     name="add_CTCF_ChIP_seq_bigWig_data_node",
            # ),
            node(
                func=find_all_motifs,
                inputs=["concat_label_HiC_loops_anotations", "readed_DNase_seq_peaks", 'motifs_JASPAR_vertebrates', 'hg19', "params:radius"],
                outputs="motifs_found",
                name="find_motifs_node",
            ),
        ]
    )
