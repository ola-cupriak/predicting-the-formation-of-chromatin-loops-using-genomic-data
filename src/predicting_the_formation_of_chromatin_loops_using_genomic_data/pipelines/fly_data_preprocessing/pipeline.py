"""
This is a boilerplate pipeline 'fly_data_preprocessing'
generated using Kedro 0.18.6
"""

from kedro.pipeline import Pipeline, node, pipeline
from .nodes import fly_read_hic, fly_read_bigWig
from .nodes import fly_add_labels
from .nodes import fly_add_bigWig_data
from .nodes import fly_all_anchors2one_df
from .nodes import fly_get_regions_with_names
from .nodes import fly_getfasta_bedfile
from .nodes import fly_find_motifs, fly_count_motifs
from .nodes import fly_remove_overlapping
from .nodes import fly_concat_dfs_from_dict


def create_pipeline(neg_sampling_type: str, **kwargs) -> Pipeline:
    namespace = neg_sampling_type

    pipeline_template = pipeline([
            node(
                func=fly_read_hic,
                inputs=["FLY_HiC_loops_annoatations", "FLY_cells2names", "params:fly_HiC_data", "params:fly_radius", "params:fly_cell_types_to_use"],
                outputs="FLY_readed_HiC_loops_anotations",
                name="FLY_read_HiC_loops_anotations_node",
            ),
            node(
                func=fly_read_bigWig,
                inputs=["FLY_DNase_seq_4-6h_bigWig", "FLY_cells2names", "params:DNase_seq_4-6h_bigWig", "params:fly_cell_types_to_use"],
                outputs="readed_DNase_seq_4-6h_bigWig",
                name="read_DNase_seq_4-6h_bigWig_data_node",
            ),
            node(
                func=fly_read_bigWig,
                inputs=["FLY_DNase_seq_6-8h_bigWig", "FLY_cells2names", "params:DNase_seq_6-8h_bigWig", "params:fly_cell_types_to_use"],
                outputs="readed_DNase_seq_6-8h_bigWig",
                name="read_DNase_seq_6-8h_bigWig_data_node",
            ),
            node(
                func=fly_read_bigWig,
                inputs=["FLY_DNase_seq_8-10h_bigWig", "FLY_cells2names", "params:DNase_seq_8-10h_bigWig", "params:fly_cell_types_to_use"],
                outputs="readed_DNase_seq_8-10h_bigWig",
                name="read_DNase_seq_8-10h_bigWig_data_node",
            ),
            node(
                func=fly_read_bigWig,
                inputs=["FLY_DNase_seq_10-12h_bigWig", "FLY_cells2names", "params:DNase_seq_10-12h_bigWig", "params:fly_cell_types_to_use"],
                outputs="readed_DNase_seq_10-12h_bigWig",
                name="read_DNase_seq_10-12h_bigWig_data_node",
            ),
            node(
                func=fly_read_bigWig,
                inputs=["FLY_Beaf-32_ChIP_seq_bigWig", "FLY_cells2names", "params:Beaf-32_ChIP_seq_bigWig", "params:fly_cell_types_to_use"],
                outputs="readed_Beaf-32_ChIP_seq_bigWig",
                name="read_Beaf-32_ChIP_seq_bigWig_data_node",
            ),
            node(
                func=fly_read_bigWig,
                inputs=["FLY_Cap-H2_ChIP_seq_bigWig", "FLY_cells2names", "params:Cap-H2_ChIP_seq_bigWig", "params:fly_cell_types_to_use"],
                outputs="readed_Cap-H2_ChIP_seq_bigWig",
                name="read_Cap-H2_ChIP_seq_bigWig_data_node",
            ),
            node(
                func=fly_read_bigWig,
                inputs=["FLY_Chromator_ChIP_seq_bigWig", "FLY_cells2names", "params:Chromator_ChIP_seq_bigWig", "params:fly_cell_types_to_use"],
                outputs="readed_Chromator_ChIP_seq_bigWig",
                name="read_Chromator_ChIP_seq_bigWig_data_node",
            ),
            node(
                func=fly_read_bigWig,
                inputs=["FLY_CP190_ChIP_seq_bigWig", "FLY_cells2names", "params:CP190_ChIP_seq_bigWig", "params:fly_cell_types_to_use"],
                outputs="readed_CP190_ChIP_seq_bigWig",
                name="read_CP190_ChIP_seq_bigWig_data_node",
            ),
            node(
                func=fly_read_bigWig,
                inputs=["FLY_CTCF_ChIP_seq_bigWig", "FLY_cells2names", "params:fly_CTCF_ChIP_seq_bigWig", "params:fly_cell_types_to_use"],
                outputs="FLY_readed_CTCF_ChIP_seq_bigWig",
                name="FLY_read_CTCF_ChIP_seq_bigWig_data_node",
            ),
            node(
                func=fly_read_bigWig,
                inputs=["FLY_GAF_ChIP_seq_bigWig", "FLY_cells2names", "params:GAF_ChIP_seq_bigWig", "params:fly_cell_types_to_use"],
                outputs="readed_GAF_ChIP_seq_bigWig",
                name="read_GAF_ChIP_seq_bigWig_data_node",
            ),
            node(
                func=fly_read_bigWig,
                inputs=["FLY_Ibf2_ChIP_seq_bigWig", "FLY_cells2names", "params:Ibf2_ChIP_seq_bigWig", "params:fly_cell_types_to_use"],
                outputs="readed_Ibf2_ChIP_seq_bigWig",
                name="read_Ibf2_ChIP_seq_bigWig_data_node",
            ),
            node(
                func=fly_read_bigWig,
                inputs=["FLY_M1BP_ChIP_seq_bigWig", "FLY_cells2names", "params:M1BP_ChIP_seq_bigWig", "params:fly_cell_types_to_use"],
                outputs="readed_M1BP_ChIP_seq_bigWig",
                name="read_M1BP_ChIP_seq_bigWig_data_node",
            ),
            node(
                func=fly_read_bigWig,
                inputs=["FLY_Pita_ChIP_seq_bigWig", "FLY_cells2names", "params:Pita_ChIP_seq_bigWig", "params:fly_cell_types_to_use"],
                outputs="readed_Pita_ChIP_seq_bigWig",
                name="read_Pita_ChIP_seq_bigWig_data_node",
            ),
            node(
                func=fly_read_bigWig,
                inputs=["FLY_Rad21_ChIP_seq_bigWig", "FLY_cells2names", "params:Rad21_ChIP_seq_bigWig", "params:fly_cell_types_to_use"],
                outputs="readed_Rad21_ChIP_seq_bigWig",
                name="read_Rad21_ChIP_seq_bigWig_data_node",
            ),
            node(
                func=fly_read_bigWig,
                inputs=["FLY_SuHw_ChIP_seq_bigWig", "FLY_cells2names", "params:SuHw_ChIP_seq_bigWig", "params:fly_cell_types_to_use"],
                outputs="readed_SuHw_ChIP_seq_bigWig",
                name="read_SuHw_ChIP_seq_bigWig_data_node",
            ),
            node(
                func=fly_read_bigWig,
                inputs=["FLY_ZIPIC_ChIP_seq_bigWig", "FLY_cells2names", "params:ZIPIC_ChIP_seq_bigWig", "params:fly_cell_types_to_use"],
                outputs="readed_ZIPIC_ChIP_seq_bigWig",
                name="read_ZIPIC_ChIP_seq_bigWig_data_node",
            ),
            node(
                func=fly_read_bigWig,
                inputs=["FLY_Zw5_ChIP_seq_bigWig", "FLY_cells2names", "params:Zw5_ChIP_seq_bigWig", "params:fly_cell_types_to_use"],
                outputs="readed_Zw5_ChIP_seq_bigWig",
                name="read_Zw5_ChIP_seq_bigWig_data_node",
            ),
            node(
                func=fly_add_labels,
                inputs=["FLY_readed_HiC_loops_anotations", "params:type", 
                        "params:fly_radius", "params:fly_neg_pos_ratio", "params:random_state"],
                outputs="FLY_positive_and_negative_HiC_loops_anotations",
                name="FLY_add_positive_and_negative_HiC_loops_anotations_node",
            ),
            node(
                func=fly_all_anchors2one_df,
                inputs=["FLY_positive_and_negative_HiC_loops_anotations"],
                outputs="FLY_merged_HiC_loops_anotations",
                name="FLY_merged_HiC_loops_anotations",
            ),
            node(
                func=fly_get_regions_with_names,
                inputs=["FLY_merged_HiC_loops_anotations"],
                outputs="FLY_overlaps_HiC_loops_DNase_seq_named",
                name="FLY_find_overlaps_HiC_loops_DNase_seq_node",
            ),
            node(
                func=fly_getfasta_bedfile,
                inputs=["FLY_overlaps_HiC_loops_DNase_seq_named", "params:path_dm6_simplified", "params:fly_path_fasta_anchors_with_open_chromtin"],
                outputs="FLY_fly_path_fasta_anchors_with_open_chromtin",
                name="FLY_getfasta_anchors_with_open_chromtin_node",
            ),
            node(
                func=fly_find_motifs,
                inputs=["params:path_motifs_JASPAR_insects", "FLY_fly_path_fasta_anchors_with_open_chromtin"],
                outputs="FLY_motifs_found_anchors_with_open_chromatin",
                name="FLY_find_motifs_in_anchors_with_open_chromtin_node",
            ),
            node(
                func=fly_add_bigWig_data,
                inputs=["FLY_positive_and_negative_HiC_loops_anotations", "readed_DNase_seq_4-6h_bigWig", "params:DNase_seq_4-6h_bigWig"],
                outputs="HiC_loops_anotations_with_DNase_seq_4-6h_bigWig_data",
                name="add_DNase_seq_4-6h_bigWig_data_node",
            ),
            node(
                func=fly_add_bigWig_data,
                inputs=["HiC_loops_anotations_with_DNase_seq_4-6h_bigWig_data", "readed_DNase_seq_6-8h_bigWig", "params:DNase_seq_6-8h_bigWig"],
                outputs="HiC_loops_anotations_with_DNase_seq_6-8h_bigWig_data",
                name="add_DNase_seq_6-8h_bigWig_data_node",
            ),
            node(
                func=fly_add_bigWig_data,
                inputs=["HiC_loops_anotations_with_DNase_seq_6-8h_bigWig_data", "readed_DNase_seq_8-10h_bigWig", "params:DNase_seq_8-10h_bigWig"],
                outputs="HiC_loops_anotations_with_DNase_seq_8-10h_bigWig_data",
                name="add_DNase_seq_8-10h_bigWig_data_node",
            ),
            node(
                func=fly_add_bigWig_data,
                inputs=["HiC_loops_anotations_with_DNase_seq_8-10h_bigWig_data", "readed_DNase_seq_10-12h_bigWig", "params:DNase_seq_10-12h_bigWig"],
                outputs="HiC_loops_anotations_with_DNase_seq_10-12h_bigWig_data",
                name="add_DNase_seq_10-12h_bigWig_data_node",
            ),
            node(
                func=fly_add_bigWig_data,
                inputs=["HiC_loops_anotations_with_DNase_seq_10-12h_bigWig_data", "readed_Beaf-32_ChIP_seq_bigWig", "params:Beaf-32_ChIP_seq_bigWig"],
                outputs="HiC_loops_anotations_with_Beaf-32_ChIP_seq_bigWig_data",
                name="add_Beaf-32_ChIP_seq_bigWig_data_node",
            ),
            node(
                func=fly_add_bigWig_data,
                inputs=["HiC_loops_anotations_with_Beaf-32_ChIP_seq_bigWig_data", "readed_Cap-H2_ChIP_seq_bigWig", "params:Cap-H2_ChIP_seq_bigWig"],
                outputs="HiC_loops_anotations_with_Cap-H2_ChIP_seq_bigWig_data",
                name="add_Cap-H2_ChIP_seq_bigWig_data_node",
            ),
            node(
                func=fly_add_bigWig_data,
                inputs=["HiC_loops_anotations_with_Cap-H2_ChIP_seq_bigWig_data", "readed_Chromator_ChIP_seq_bigWig", "params:Chromator_ChIP_seq_bigWig"],
                outputs="HiC_loops_anotations_with_Chromator_ChIP_seq_bigWig_data",
                name="add_Chromator_ChIP_seq_bigWig_data_node",
            ),
            node(
                func=fly_add_bigWig_data,
                inputs=["HiC_loops_anotations_with_Chromator_ChIP_seq_bigWig_data", "readed_CP190_ChIP_seq_bigWig", "params:CP190_ChIP_seq_bigWig"],
                outputs="HiC_loops_anotations_with_CP190_ChIP_seq_bigWig_data",
                name="add_CP190_ChIP_seq_bigWig_data_node",
            ),
            node(
                func=fly_add_bigWig_data,
                inputs=["HiC_loops_anotations_with_CP190_ChIP_seq_bigWig_data", "FLY_readed_CTCF_ChIP_seq_bigWig", "params:fly_CTCF_ChIP_seq_bigWig"],
                outputs="HiC_loops_anotations_with_FLY_CTCF_ChIP_seq_bigWig_data",
                name="add_FLY_CTCF_ChIP_seq_bigWig_data_node",
            ),
            node(
                func=fly_add_bigWig_data,
                inputs=["HiC_loops_anotations_with_FLY_CTCF_ChIP_seq_bigWig_data", "readed_GAF_ChIP_seq_bigWig", "params:GAF_ChIP_seq_bigWig"],
                outputs="HiC_loops_anotations_with_GAF_ChIP_seq_bigWig_data",
                name="add_GAF_ChIP_seq_bigWig_data_node",
            ),
            node(
                func=fly_add_bigWig_data,
                inputs=["HiC_loops_anotations_with_GAF_ChIP_seq_bigWig_data", "readed_Ibf2_ChIP_seq_bigWig", "params:Ibf2_ChIP_seq_bigWig"],
                outputs="HiC_loops_anotations_with_Ibf2_ChIP_seq_bigWig_data",
                name="add_Ibf2_ChIP_seq_bigWig_data_node",
            ),
            node(
                func=fly_add_bigWig_data,
                inputs=["HiC_loops_anotations_with_Ibf2_ChIP_seq_bigWig_data", "readed_M1BP_ChIP_seq_bigWig", "params:M1BP_ChIP_seq_bigWig"],
                outputs="HiC_loops_anotations_with_M1BP_ChIP_seq_bigWig_data",
                name="add_M1BP_ChIP_seq_bigWig_data_node",
            ),
            node(
                func=fly_add_bigWig_data,
                inputs=["HiC_loops_anotations_with_M1BP_ChIP_seq_bigWig_data", "readed_Pita_ChIP_seq_bigWig", "params:Pita_ChIP_seq_bigWig"],
                outputs="HiC_loops_anotations_with_Pita_ChIP_seq_bigWig_data",
                name="add_Pita_ChIP_seq_bigWig_data_node",
            ),
            node(
                func=fly_add_bigWig_data,
                inputs=["HiC_loops_anotations_with_Pita_ChIP_seq_bigWig_data", "readed_Rad21_ChIP_seq_bigWig", "params:Rad21_ChIP_seq_bigWig"],
                outputs="HiC_loops_anotations_with_Rad21_ChIP_seq_bigWig_data",
                name="add_Rad21_ChIP_seq_bigWig_data_node",
            ),
            node(
                func=fly_add_bigWig_data,
                inputs=["HiC_loops_anotations_with_Rad21_ChIP_seq_bigWig_data", "readed_SuHw_ChIP_seq_bigWig", "params:SuHw_ChIP_seq_bigWig"],
                outputs="HiC_loops_anotations_with_SuHw_ChIP_seq_bigWig_data",
                name="add_SuHw_ChIP_seq_bigWig_data_node",
            ),
            node(
                func=fly_add_bigWig_data,
                inputs=["HiC_loops_anotations_with_SuHw_ChIP_seq_bigWig_data", "readed_ZIPIC_ChIP_seq_bigWig", "params:ZIPIC_ChIP_seq_bigWig"],
                outputs="HiC_loops_anotations_with_ZIPIC_ChIP_seq_bigWig_data",
                name="add_ZIPIC_ChIP_seq_bigWig_data_node",
            ),
            node(
                func=fly_add_bigWig_data,
                inputs=["HiC_loops_anotations_with_ZIPIC_ChIP_seq_bigWig_data", "readed_Zw5_ChIP_seq_bigWig", "params:Zw5_ChIP_seq_bigWig"],
                outputs="HiC_loops_anotations_with_Zw5_ChIP_seq_bigWig_data",
                name="add_Zw5_ChIP_seq_bigWig_data_node",
            ),
            node(
                func=fly_count_motifs,
                inputs=["HiC_loops_anotations_with_Zw5_ChIP_seq_bigWig_data", "FLY_motifs_found_anchors_with_open_chromatin"],
                outputs="FLY_combined_functional_genomics_data",
                name="FLY_add_motif_counts_node",
            ),
            node(
                func=fly_remove_overlapping,
                inputs="FLY_combined_functional_genomics_data",
                outputs="FLY_combined_functional_genomics_data_no_overlapping",
                name="FLY_remove_overlapping_node",
            ),
            node(
                func=fly_concat_dfs_from_dict,
                inputs=["FLY_combined_functional_genomics_data_no_overlapping", "params:fly_cell_types_to_use"],
                outputs="FLY_concatenated_combined_functional_genomics_data",
                name="FLY_concatenate_combined_functional_genomics_data_node",
            ),
        ])
    
    main_pipeline = pipeline(
        pipe=pipeline_template,
        inputs=[
            "FLY_HiC_loops_annoatations", 
            "FLY_cells2names", 
            "FLY_DNase_seq_4-6h_bigWig", 
            "FLY_DNase_seq_6-8h_bigWig", 
            "FLY_DNase_seq_8-10h_bigWig", 
            "FLY_DNase_seq_10-12h_bigWig", 
            "FLY_Beaf-32_ChIP_seq_bigWig",
            "FLY_Cap-H2_ChIP_seq_bigWig",
            "FLY_Chromator_ChIP_seq_bigWig",
            "FLY_CP190_ChIP_seq_bigWig",
            "FLY_CTCF_ChIP_seq_bigWig",
            "FLY_GAF_ChIP_seq_bigWig",
            "FLY_Ibf2_ChIP_seq_bigWig",
            "FLY_M1BP_ChIP_seq_bigWig",
            "FLY_Pita_ChIP_seq_bigWig",
            "FLY_Rad21_ChIP_seq_bigWig",
            "FLY_SuHw_ChIP_seq_bigWig",
            "FLY_ZIPIC_ChIP_seq_bigWig",
            "FLY_Zw5_ChIP_seq_bigWig", 
        ],
        parameters=[
            "params:fly_HiC_data",
            "params:fly_radius", 
            "params:random_state",
            "params:fly_cell_types_to_use",
            "params:DNase_seq_4-6h_bigWig",
            "params:DNase_seq_6-8h_bigWig",
            "params:DNase_seq_8-10h_bigWig",
            "params:DNase_seq_10-12h_bigWig",
            "params:Beaf-32_ChIP_seq_bigWig",
            "params:Cap-H2_ChIP_seq_bigWig",
            "params:Chromator_ChIP_seq_bigWig",
            "params:CP190_ChIP_seq_bigWig",
            "params:fly_CTCF_ChIP_seq_bigWig",
            "params:GAF_ChIP_seq_bigWig",
            "params:Ibf2_ChIP_seq_bigWig",
            "params:M1BP_ChIP_seq_bigWig",
            "params:Pita_ChIP_seq_bigWig",
            "params:Rad21_ChIP_seq_bigWig",
            "params:SuHw_ChIP_seq_bigWig",
            "params:ZIPIC_ChIP_seq_bigWig",
            "params:Zw5_ChIP_seq_bigWig",
            "params:path_dm6_simplified", 
            "params:path_motifs_JASPAR_insects",
        ],
        namespace=namespace,
    )

    return main_pipeline

