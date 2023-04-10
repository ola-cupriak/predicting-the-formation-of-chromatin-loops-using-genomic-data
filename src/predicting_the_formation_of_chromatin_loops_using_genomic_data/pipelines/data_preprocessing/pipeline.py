"""
This is a boilerplate pipeline 'data_preprocessing'
generated using Kedro 0.18.6
"""

from kedro.pipeline import Pipeline, node, pipeline
from .nodes import read_hic, read_peaks, read_bigWig
from .nodes import add_labels
from .nodes import count_peaks
from .nodes import add_bigWig_data
from .nodes import all_anchors2one_df, all_peaks2one_df
from .nodes import get_overlapping_regions, get_overlaps_with_names
from .nodes import getfasta_bedfile
from .nodes import find_motifs, count_motifs


def create_pipeline(**kwargs) -> Pipeline:
    return pipeline(
        [
            node(
                func=read_hic,
                inputs=["HiC_loops_annoatations", "cells2names", "params:HiC_data", "params:radius"],
                outputs="readed_HiC_loops_anotations",
                name="read_HiC_loops_anotations_node",
            ),
            # node(
            #     func=read_peaks,
            #     inputs=["DNAse_seq_peaks", "cells2names", "params:DNase-seq_peaks"],
            #     outputs="readed_DNase_seq_peaks",
            #     name="read_DNase_seq_peaks_node",
            # ),
            node(
                func=add_labels,
                inputs="readed_HiC_loops_anotations",
                outputs="concat_label_HiC_loops_anotations",
                name="concat_label_HiC_loops_anotations_node",
            ),
            # node(
            #     func=all_anchors2one_df,
            #     inputs=["readed_HiC_loops_anotations", "params:radius"],
            #     outputs="merged_HiC_loops_anotations",
            #     name="merge_HiC_loops_anotations_to_one_df_node",
            # ),
            # node(
            #     func=all_peaks2one_df,
            #     inputs=["readed_DNase_seq_peaks"],
            #     outputs="merged_DNase_seq_peaks",
            #     name="merge_DNase_seq_peaks_to_one_df_node",
            # ),
            # node(
            #     func=get_overlaps_with_names,
            #     inputs=["merged_HiC_loops_anotations", "merged_DNase_seq_peaks"],
            #     outputs="overlaps_HiC_loops_DNase_seq_named",
            #     name="find_overlaps_HiC_loops_DNase_seq_node",
            # ),
            # node(
            #     func=getfasta_bedfile,
            #     inputs=["overlaps_HiC_loops_DNase_seq_named", "params:path_hg19_simplified", "params:path_fasta_anchors_with_open_chromtin"],
            #     outputs="path_fasta_anchors_with_open_chromtin",
            #     name="getfasta_anchors_with_open_chromtin_node",
            # ),
            # node(
            #     func=find_motifs,
            #     inputs=["params:path_motifs_JASPAR_vertebrates", "params:path_fasta_anchors_with_open_chromtin"],
            #     outputs="motifs_found_anchors_with_open_chromatin",
            #     name="find_motifs_in_anchors_with_open_chromtin_node",
            # ),
            node(
                func=count_motifs,
                inputs=["concat_label_HiC_loops_anotations", "motifs_found_anchors_with_open_chromatin"],
                outputs="combined_functional_genomics_data",
                name="add_motif_counts_node",
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
        ]
    )
