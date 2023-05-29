import pandas as pd
import polars as pl
from predicting_the_formation_of_chromatin_loops_using_genomic_data.utils import _dict_partitions
from predicting_the_formation_of_chromatin_loops_using_genomic_data.pipelines.data_preprocessing.nodes import add_labels
from predicting_the_formation_of_chromatin_loops_using_genomic_data.pipelines.data_preprocessing.nodes import read_bigWig, add_bigWig_data
from predicting_the_formation_of_chromatin_loops_using_genomic_data.pipelines.data_preprocessing.nodes import all_anchors2one_df, getfasta_bedfile
from predicting_the_formation_of_chromatin_loops_using_genomic_data.pipelines.data_preprocessing.nodes import find_motifs, count_motifs
from predicting_the_formation_of_chromatin_loops_using_genomic_data.pipelines.data_preprocessing.nodes import remove_overlapping, concat_dfs_from_dict
from typing import Any, Callable, Dict



def _prepare_fly_loops_data(df: pd.DataFrame, r: int) -> pd.DataFrame:
    """
    Create dataframe with pairs of anchors creating chromatin loops for fly.
    Args:
        df: pandas DataFrame with chromatin loops anotations (each line is a single anchor).
        r: radius of the region.
    Returns:
        pandas DataFrame with pairs of anchors creating chromatin loops.
    """
    df = df[['loop_id', 'anchor', 'anchor_chr', 'anchor_summit']]
    df = pl.from_pandas(df)
    df = df.join(df, on='loop_id', suffix='_y').filter(pl.col('anchor_summit') != pl.col('anchor_summit_y')).filter(pl.col('anchor') == 1).filter(pl.col('anchor_y') == 2).filter(pl.col('anchor_chr') == pl.col('anchor_chr_y'))
    df = df.to_pandas()
    df = df[['anchor_chr', 'anchor_summit', 'anchor_summit_y']]
    df.rename(columns={'anchor_chr': 'chr', 'anchor_summit': 'x', 'anchor_summit_y': 'y'}, inplace=True)
    df['x_start'] = df['x'] - r
    df['x_end'] = df['x'] + r
    df['y_start'] = df['y'] - r
    df['y_end'] = df['y'] + r
    df = df.astype({'x': 'int32', 'y': 'int32', 'x_start': 'int32', 'x_end': 'int32', 
                    'y_start': 'int32', 'y_end': 'int32', 'chr': 'string'})
    df = df = df.sort_values(by=['chr', 'x'])
    
    return df

def fly_read_hic(partitioned_input: Dict[str, Callable[[], Any]], 
                cells2names: Dict[str, dict],
                dataset_name: str, r: int,
                cells_to_use: list=[]) -> Dict[str, pd.DataFrame]:
    """
    Load and modify the dataframes with chromatin loops anotations for fly.
    Args:
        partitioned_input: dictionary with partition ids as keys and load functions as values.
        cells2names: dictionary, template: {'dataset_name': {'file_name': 'cell_type'}}
        dataset_name: name of the dataset to select from cells2names.
        r: radius of the region.
        cells_to_use: list of cell types to use.
    Returns:
        dictionary:
            keys: cell types 
            values: pandas DataFrames with chromatin loops anotations.
    """
    dfs_dict = _dict_partitions(partitioned_input)
    cells2names_dataset_dict = cells2names[dataset_name]
    assert set(cells_to_use).issubset(set(cells2names_dataset_dict.values())), f"Cell types: {set(cells_to_use)-set(cells2names_dataset_dict.values())} are not in the dataset. Please check data_preprocessing.yml file."
    new_dfs_dict = dict()
    for name, df in dfs_dict.items():
        cell_type = cells2names_dataset_dict[name]
        if cells_to_use and cell_type not in cells_to_use:
            continue
        df = _prepare_fly_loops_data(df, r)
        # set dtypes
        df['cell_type'] = cell_type
        new_dfs_dict[cell_type] = df

    return new_dfs_dict


def fly_add_labels(dfs_dict: Dict[str, pd.DataFrame], mtype: str, r, neg_pos_ratio: float, random_state: int, peaks_dict: Dict[str, pd.DataFrame]=None) -> None:
    """
    Add labels to the dataframes depending on the cell type and type of model to train.
    Within model: 1 if cell type is the same as the cell type of the loop, 0 otherwise.
    Across model: negative sampling involving the use of open regions of chromatin that 
                do not overlap with any positive example
    Args:
        dfs_dict: dictionary with cell types as keys and pandas DataFrames as values.
        type: type of model to train, either 'within' or 'across'.
        peaks_dict: dictionary with cell types as keys and pandas DataFrames with peaks as values.
        r
    Returns:
        dictionary:
            keys: cell types
            values: pandas DataFrames with labels added.
    """
    return add_labels(dfs_dict, mtype, peaks_dict, r, neg_pos_ratio, random_state)


def fly_read_bigWig(partitioned_input: Dict[str, Callable[[], Any]],
                cells2names: Dict[str, dict],
                dataset_name: str,
                cells_to_use: list) -> Dict[str, pd.DataFrame]:
    """
    Create a dictionary with paths to DNase-seq/CTCF ChIP-seq bigWig files for each cell type.
    Args:
        partitioned_input: dictionary with partition ids as keys and load functions as values.
        cells2names: dictionary, template: {'dataset_name': {'file_name': 'cell_type'}}
        dataset_name: name of the dataset to select from cells2names.
        cells_to_use: list of cell types to use.
    Returns:
        dictionary:
            keys: cell types 
            values: paths to DNase-seq/CTCF ChIP-seq bigWig files.
    """
    return read_bigWig(partitioned_input, cells2names, dataset_name, cells_to_use)


def fly_add_bigWig_data(main_dfs_dict: Dict[str, Callable[[], Any]],
                    bigWig_data_dict: dict,
                    experiment: str) -> pd.DataFrame:
    """
    Count statistics (weighted mean, arithmetic mean, minimum and maximum) 
    of the bigWig data in both regions of each chromatin loop,
    for each dataframe from the main_dfs_dict dictionary.
    Args:
        main_dfs_dict: dictionary with cell types as keys and load functions of pandas DataFrames with chromatin loops as values.
        bigWig_data_dict: dictionary with cell types as keys and pyBigWig objects as values.
    Returns:
        dictionary:
            keys: cell types 
            values: pandas DataFrames with added columns of bigWig data statistics in both regions of each loop
    """
    return add_bigWig_data(main_dfs_dict, bigWig_data_dict, experiment, organism='fly')


def fly_all_anchors2one_df(dfs_dict: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    """
    Concatenates all anchors DataFrames from dfs_dict into one DataFrame and 
    combines the columns describing x regions and the columns describing y regions 
    into one set of columns describing all regions.
    (columns: x_chr, x_start, x_end, y_chr, y_start, y_end -> columns: chr, start, end)
    Args:
        dfs_dict: dictionary with cell types as keys and pandas DataFrames as values.
    Returns:
        pandas DataFrame with one set of columns describing all regions.
    """
    return all_anchors2one_df(dfs_dict)


def fly_get_regions_with_names(anchors_df: pd.DataFrame) -> pd.DataFrame:
    """
    Get overlapping regions between anchors and peaks.
    Args:
        anchors_df: pandas DataFrame with anchors coordinates.
    Returns:
        pandas DataFrame with overlapping regions coordinates and names.
    """
    anchors_df['name'] = anchors_df.apply(lambda x: f"{x['chr']}:{x['start']}:{x['end']}:neuronal", axis=1)

    return anchors_df[['chr', 'start', 'end', 'name']]


def fly_getfasta_bedfile(df: pd.DataFrame, path_simp_genome: str, path_to_save: str) -> str:
    """
    Cut sequences from chromosomes using BEDTools for coordinates from the pandas DataFrame.
    Args:
        df: pandas DataFrame with coordinates.
        path_simp_genome: path to fasta file with chromosomes.
        path_to_save: path to save fasta file.
    Returns:
        path_to_save: path to fasta file.
    """
    return getfasta_bedfile(df, path_simp_genome, path_to_save, organism='fly')


def fly_find_motifs(path_motifs: str, path_fasta: list) -> pd.DataFrame:
    """
    Finds motif occurances in fasta sequences.
    Args:
        path_motifs: path to file with motifs in jaspar format.
        path_fasta: path to fasta file with sequences.
    Returns:
        pandas DataFrame with counts of each motif found in fasta sequences.
    """
    return find_motifs(path_motifs, path_fasta)


def fly_count_motifs(main_dfs_dict: dict, motifs_df: pd.DataFrame):
    """
    Counts motif occurances in both regions of each chromatin loop,
    for each dataframe from the main_dfs_dict dictionary.
    Args:
        main_dfs_dict: dictionary with cell types as keys and load functions of pandas DataFrames with chromatin loops as values.
        motifs_df: dictionary with cell types as keys and pandas DataFrame with motif occurrences found as values.
    Returns:
        dictionary:
            keys: cell types 
            values: pandas DataFrame with added columns of each motif counts in both regions of each chromatin loop
    """
    return count_motifs(main_dfs_dict, motifs_df)


def fly_remove_overlapping(main_dfs_dict: dict):
    """
    Remove negative examples that overlap with positive examples,
    for each dataframe from the main_dfs_dict dictionary.
    Args:
        main_dfs_dict: dictionary with cell types as keys and load functions of pandas DataFrames with chromatin loops as values.
    Returns:
        dictionary:
            keys: cell types
            values: pandas DataFrame with removed negative examples that overlap with positive examples
    """
    return remove_overlapping(main_dfs_dict)


def fly_concat_dfs_from_dict(main_dfs_dict: dict, cells_to_use: list=[]) -> pd.DataFrame:
    '''
    Concatenates dataframes from dictionary.
    Args:   
        main_dfs_dict: dictionary with cell types as keys and load functions of pandas DataFrames with chromatin loops as values.
        cells_to_use: list of cell types to be used. If empty, all cell types from main_dfs_dict will be used.
    Returns:
        pandas DataFrame with concatenated dataframes from dictionary.
    '''
    return concat_dfs_from_dict(main_dfs_dict, cells_to_use)