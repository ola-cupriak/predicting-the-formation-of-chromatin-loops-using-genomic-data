"""
This is a boilerplate pipeline 'data_preprocessing'
generated using Kedro 0.18.6
"""

from kedro.pipeline import Pipeline, node, pipeline
from .nodes import read_hic, read_peaks, read_bigWig
from .nodes import add_labels
from .nodes import count_peaks_and_distances
from .nodes import add_bigWig_data
from .nodes import all_anchors2one_df, all_peaks2one_df
from .nodes import get_overlaps_with_names
from .nodes import getfasta_bedfile
from .nodes import find_motifs, count_motifs
from .nodes import remove_overlapping
from .nodes import concat_dfs_from_dict


def create_pipeline(neg_sampling_type: str, **kwargs) -> Pipeline:
    namespace = neg_sampling_type

    pipeline_template =  pipeline([
            node(
                func=read_hic,
                inputs=["HiC_loops_annoatations", "cells2names", "params:HiC_data", "params:radius", "params:cell_types_to_use"],
                outputs="readed_HiC_loops_anotations",
                name="read_HiC_loops_anotations_node",
            ),
            node(
                func=read_peaks,
                inputs=["DNAse_seq_peaks", "cells2names", "params:DNase-seq_peaks", "params:cell_types_to_use"],
                outputs="readed_DNase_seq_peaks",
                name="read_DNase_seq_peaks_node",
            ),
            node(
                func=read_peaks,
                inputs=["CTCF_ChIP_seq_peaks", "cells2names", "params:CTCF_ChIP-seq_peaks", "params:cell_types_to_use"],
                outputs="readed_CTCF_ChIP_seq_peaks",
                name="read_CTCF_ChIP_seq_peaks_node",
            ),
            node(
                func=read_bigWig,
                inputs=["DNAse_seq_bigWig", "cells2names", "params:DNase-seq_bigWig", "params:cell_types_to_use"],
                outputs="readed_DNase_seq_bigWig",
                name="read_DNase_seq_bigWig_data_node",
            ),
            node(
                func=read_bigWig,
                inputs=["CTCF_ChIP_seq_bigWig", "cells2names", "params:CTCF_ChIP-seq_bigWig", "params:cell_types_to_use"],
                outputs="readed_CTCF_ChIP_seq_bigWig",
                name="read_CTCF_ChIP_seq_bigWig_data_node",
            ),
            node(
                func=add_labels,
                inputs=["readed_HiC_loops_anotations", "params:type", "readed_DNase_seq_peaks", 
                        "params:radius", "params:human_neg_pos_ratio", "params:random_state"],
                outputs="positive_and_negative_HiC_loops_anotations",
                name="add_positive_and_negative_HiC_loops_anotations_node",
            ),
            node(
                func=all_anchors2one_df,
                inputs=["positive_and_negative_HiC_loops_anotations"],
                outputs="merged_HiC_loops_anotations",
                name="merge_HiC_loops_anotations_to_one_df_node",
            ),
            node(
                func=all_peaks2one_df,
                inputs=["readed_DNase_seq_peaks"],
                outputs="merged_DNase_seq_peaks",
                name="merge_DNase_seq_peaks_to_one_df_node",
            ),
            node(
                func=get_overlaps_with_names,
                inputs=["merged_HiC_loops_anotations", "merged_DNase_seq_peaks"],
                outputs="overlaps_HiC_loops_DNase_seq_named",
                name="find_overlaps_HiC_loops_DNase_seq_node",
            ),
            node(
                func=getfasta_bedfile,
                inputs=["overlaps_HiC_loops_DNase_seq_named", "params:path_hg19_simplified", "params:human_path_fasta_anchors_with_open_chromtin"],
                outputs="human_path_fasta_anchors_with_open_chromtin",
                name="getfasta_anchors_with_open_chromtin_node",
            ),
            node(
                func=find_motifs,
                inputs=["params:path_motifs_JASPAR_vertebrates", "human_path_fasta_anchors_with_open_chromtin"],
                outputs="motifs_found_anchors_with_open_chromatin",
                name="find_motifs_in_anchors_with_open_chromtin_node",
            ),
            node(
                func=count_peaks_and_distances,
                inputs=["positive_and_negative_HiC_loops_anotations", "readed_CTCF_ChIP_seq_peaks", "params:CTCF_ChIP-seq_peaks"],
                outputs="HiC_loops_anotations_with_CTCF_ChIP_seq_peaks",
                name="add_CTCF_ChIP_seq_peaks_node",
            ),
            node(
                func=count_peaks_and_distances,
                inputs=["HiC_loops_anotations_with_CTCF_ChIP_seq_peaks", "readed_DNase_seq_peaks", "params:DNase-seq_peaks"],
                outputs="HiC_loops_anotations_with_DNase_seq_peaks",
                name="add_DNase_seq_peaks_node",
            ),
            node(
                func=add_bigWig_data,
                inputs=["HiC_loops_anotations_with_DNase_seq_peaks", "readed_DNase_seq_bigWig", "params:DNase-seq_bigWig", "params:resolution"],
                outputs="HiC_loops_anotations_with_DNase_seq_bigWig_data",
                name="add_DNase_seq_bigWig_data_node",
            ),
            node(
                func=add_bigWig_data,
                inputs=["HiC_loops_anotations_with_DNase_seq_bigWig_data", "readed_CTCF_ChIP_seq_bigWig", "params:CTCF_ChIP-seq_bigWig", "params:resolution"],
                outputs="HiC_loops_anotations_with_CTCF_ChIP_seq_bigWig_data",
                name="add_CTCF_ChIP_seq_bigWig_data_node",
            ),
            node(
                func=count_motifs,
                inputs=["HiC_loops_anotations_with_CTCF_ChIP_seq_bigWig_data", "motifs_found_anchors_with_open_chromatin"],
                outputs="combined_functional_genomics_data",
                name="add_motif_counts_node",
            ),
            node(
                func=remove_overlapping,
                inputs="combined_functional_genomics_data",
                outputs="combined_functional_genomics_data_no_overlapping",
                name="remove_overlapping_node",
            ),
            node(
                func=concat_dfs_from_dict,
                inputs=["combined_functional_genomics_data_no_overlapping", "params:cell_types_to_use"],
                outputs="concatenated_combined_functional_genomics_data",
                name="concatenate_combined_functional_genomics_data_node",
            ),
        ])
    
    main_pipeline = pipeline(
        pipe=pipeline_template,
        inputs=[
            "HiC_loops_annoatations", 
            "cells2names", 
            "DNAse_seq_peaks", 
            "CTCF_ChIP_seq_peaks", 
            "DNAse_seq_bigWig", 
            "CTCF_ChIP_seq_bigWig"
        ],
        parameters=[
            "params:HiC_data",
            "params:radius", 
            "params:resolution",
            "params:random_state",
            "params:cell_types_to_use",
            "params:DNase-seq_peaks",
            "params:CTCF_ChIP-seq_peaks",
            "params:DNase-seq_bigWig",
            "params:CTCF_ChIP-seq_bigWig",
            "params:path_hg19_simplified", 
            "params:path_motifs_JASPAR_vertebrates",
        ],
        namespace=namespace,
    )

    return main_pipeline

