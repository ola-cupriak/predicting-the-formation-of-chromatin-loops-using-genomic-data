"""
This is a boilerplate pipeline 'data_preprocessing_old'
generated using Kedro 0.18.6
"""

from kedro.pipeline import Pipeline, node, pipeline
from .nodes import read_hic, add_labels, read_peaks, count_peaks, read_bigWig, add_bigWig_data, gather_all_anchors_into_df, getfasta_anchors, anchors2df

def create_pipeline(**kwargs) -> Pipeline:
    return pipeline(
        [
            node(
                func=read_hic,
                inputs=["HiC_loops_annoatations", "cells2names", "params:HiC_data", 'params:radius'],
                outputs="readed_HiC_loops_anotations",
                name="read_HiC_loops_anotations_node",
            ),
            # node(
            #     func=gather_all_anchors_into_df,
            #     inputs=["readed_HiC_loops_anotations", "params:radius"],
            #     outputs=["gathered_HiC_loops_anotations", "path_gathered_HiC_loops_anotations"],
            #     name="gather_HiC_loops_anotations_into_df_node",
            # ),
            # node(
            #     func=getfasta_anchors,
            #     inputs=["path_gathered_HiC_loops_anotations", "path_hg19_simplified"],
            #     outputs="fasta_anchors",
            #     name="getfasta_anchors_node",
            # ),
            # node(
            #     func=anchors2df,
            #     inputs="fasta_anchors",
            #     outputs="anchors_df",
            #     name="fasta2df_anchor_node",
            # ),
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
            node(
                func=read_peaks,
                inputs=["CTCF_ChIP_seq_peaks", "cells2names", "params:CTCF_ChIP-seq_peaks"],
                outputs="readed_CTCF_ChIP_seq_peaks",
                name="read_CTCF_ChIP_seq_peaks_node",
            ),
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
            node(
                func=count_peaks,
                inputs=["concat_label_HiC_loops_anotations", "readed_CTCF_ChIP_seq_peaks", "anchors_df", "params:CTCF_ChIP-seq_peaks", "params:radius"],
                outputs="HiC_loops_anotations_with_CTCF_ChIP_seq_peaks",
                name="add_CTCF_ChIP_seq_peaks_node",
            ),
            node(
                func=count_peaks,
                inputs=["HiC_loops_anotations_with_CTCF_ChIP_seq_peaks", "readed_DNase_seq_peaks", "anchors_df", "params:DNase-seq_peaks", "params:radius"],
                outputs="combined_functional_genomics_data",
                name="add_DNase_seq_peaks_and_motifs_node",
            ),
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
        ]
    )
