import pandas as pd
import polars as pl
from predicting_the_formation_of_chromatin_loops_using_genomic_data.utils import (
    _dict_partitions,
)
from predicting_the_formation_of_chromatin_loops_using_genomic_data.pipelines.data_preprocessing.nodes import (
    add_labels,
    _remove_overlapping_single_df,
)
from predicting_the_formation_of_chromatin_loops_using_genomic_data.pipelines.data_preprocessing.nodes import (
    read_peaks,
    count_peaks_and_distances,
)
from predicting_the_formation_of_chromatin_loops_using_genomic_data.pipelines.data_preprocessing.nodes import (
    read_bigWig,
    add_bigWig_data,
)
from predicting_the_formation_of_chromatin_loops_using_genomic_data.pipelines.data_preprocessing.nodes import (
    all_anchors2one_df,
    getfasta_bedfile,
    all_peaks2one_df,
    get_overlaps_with_names,
)
from predicting_the_formation_of_chromatin_loops_using_genomic_data.pipelines.data_preprocessing.nodes import (
    find_motifs,
    count_motifs,
)
from predicting_the_formation_of_chromatin_loops_using_genomic_data.pipelines.data_preprocessing.nodes import (
    remove_overlapping,
    concat_dfs_from_dict,
)
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
    df = df[["chr1", "x1", "x2", "y1", "y2"]]
    df["x"] = (df["x1"] + df["x2"]) // 2
    df["y"] = (df["y1"] + df["y2"]) // 2
    df.rename(
        columns={"chr1": "chr"},
        inplace=True,
    )
    df["x_start"] = df["x"] - r
    df["x_end"] = df["x"] + r
    df["y_start"] = df["y"] - r
    df["y_end"] = df["y"] + r
    df.drop(["x1", "x2", "y1", "y2"], axis=1, inplace=True)
    df = df.astype(
        {
            "x": "int32",
            "y": "int32",
            "x_start": "int32",
            "x_end": "int32",
            "y_start": "int32",
            "y_end": "int32",
            "chr": "string",
        }
    )
    df = df = df.sort_values(by=["chr", "x"])

    return df


def fly_read_hic(
    partitioned_input: Dict[str, Callable[[], Any]],
    cells2names: Dict[str, dict],
    dataset_name: str,
    r: int,
) -> Dict[str, pd.DataFrame]:
    """
    Load and modify the dataframes with chromatin loops anotations for fly.
    Args:
        partitioned_input: dictionary with partition ids as keys and load functions as values.
        cells2names: dictionary, template: {'dataset_name': {'file_name': 'cell_type'}}
        dataset_name: name of the dataset to select from cells2names.
        r: radius of the region.
    Returns:
        dictionary:
            keys: cell types
            values: pandas DataFrames with chromatin loops anotations.
    """
    dfs_dict = _dict_partitions(partitioned_input)
    cells2names_dataset_dict = cells2names[dataset_name]
    new_dfs_dict = dict()
    for name, df in dfs_dict.items():
        cell_type = cells2names_dataset_dict[name]
        df = _prepare_fly_loops_data(df, r)
        # set dtypes
        df["cell_type"] = cell_type
        new_dfs_dict[cell_type] = df

    return new_dfs_dict


def fly_read_peaks(
    partitioned_input: Dict[str, Callable[[], Any]],
    cells2names: Dict[str, dict],
    dataset_name: str,
) -> Dict[str, pd.DataFrame]:
    """
    Load dataframes with experiment peaks saved in bed file and modify the chromosome columns.
    Args:
        partitioned_input: dictionary with partition ids as keys and load functions as values.
        cells2names: dictionary, template: {'dataset_name': {'file_name': 'cell_type'}}
        dataset_name: name of the dataset to select from cells2names.
    Returns:
        dictionary:
            keys: cell types
            values: pandas DataFrames with DNase-seq/ChIP-seq peaks.
    """
    return read_peaks(
        partitioned_input, cells2names, dataset_name, cells_to_use=[], organism="fly"
    )


def fly_read_bigWig(
    partitioned_input: Dict[str, Callable[[], Any]],
    cells2names: Dict[str, dict],
    dataset_name: str,
) -> Dict[str, pd.DataFrame]:
    """
    Create a dictionary with paths to DNase-seq/CTCF ChIP-seq bigWig files for each cell type.
    Args:
        partitioned_input: dictionary with partition ids as keys and load functions as values.
        cells2names: dictionary, template: {'dataset_name': {'file_name': 'cell_type'}}
        dataset_name: name of the dataset to select from cells2names.
    Returns:
        dictionary:
            keys: cell types
            values: paths to DNase-seq/CTCF ChIP-seq bigWig files.
    """
    return read_bigWig(partitioned_input, cells2names, dataset_name, cells_to_use=[])


def _neg_sampling_short_loops(df: pd.DataFrame, short_long_limit: int) -> pd.DataFrame:
    """ """
    df["label"] = [1 if diff > short_long_limit else 0 for diff in df["y"] - df["x"]]
    len_before = len(df)
    df = _remove_overlapping_single_df(df)
    print(f"Removed {len_before - len(df)} overlapping negative loops.")
    return df


def fly_add_labels(
    dfs_dict: Dict[str, pd.DataFrame],
    negtype: str,
    peaks_dict: Dict[str, pd.DataFrame],
    r: int,
    neg_pos_ratio: float,
    random_state: int,
    short_long_limit: int = 1000000,
) -> Dict[str, pd.DataFrame]:
    """
    Add labels to the dataframes depending on the cell type and type of model to train.
    Within model: 1 if cell type is the same as the cell type of the loop, 0 otherwise.
    Across model: negative sampling involving the use of open regions of chromatin that
                do not overlap with any positive example
    Args:
        dfs_dict: dictionary with cell types as keys and pandas DataFrames as values.
        negtype: type of negative sampling to use.
        peaks_dict: dictionary with cell types as keys and pandas DataFrames with peaks as values.
        r: radius of the region around the anchor.
        neg_pos_ratio: ratio of negative to positive examples.
        random_state: random state.
        short_long_limit: limit for short loops.
    Returns:
        dictionary:
            keys: cell types
            values: pandas DataFrames with labels added.
    """
    if negtype in [
        "anchors_from_other_cell_types",
        "new_pairs_of_anchors",
        "open_chromatin_regions",
    ]:
        return add_labels(
            dfs_dict,
            negtype,
            peaks_dict,
            r,
            neg_pos_ratio,
            random_state,
            organism="fly",
        )
    elif negtype == "FLY_short_loops_as_negatives":
        for name, cell_df in dfs_dict.items():
            print("Creating negative examples for cell type", name, "...")
            cell_df = _neg_sampling_short_loops(cell_df, short_long_limit)
            dfs_dict[name] = cell_df
        return dfs_dict
    else:
        raise ValueError("Wrong negative sampling type.")


def fly_count_peaks_and_distances(
    main_dfs_dict: Dict[str, Callable[[], Any]],
    peaks_dfs_dict: Dict[str, pd.DataFrame],
    experiment: str,
) -> pd.DataFrame:
    """
    Count the number of peaks in both regions of each chromatin loop
    and find the distance of the closest peak from the center of each anchor from the chromatin loop,
    for each dataframe from the main_dfs_dict dictionary.
    Args:
        main_dfs_dict: dictionary with cell types as keys and load functions of pandas DataFrames with chromatin loops as values.
        peaks_dfs_dict: dictionary with cell types as keys and pandas DataFrames with experiment peaks as values.
        experiment: name of the experiment.
    Returns:
        dictionary:
            keys: cell types
            values: pandas DataFrame with added columns of numbers of experiment peaks in both regions of each chromatin loop
                    and columns of distances of the closest peak of the experiment from the center of each anchor from the chromatin loop.
    """
    return count_peaks_and_distances(
        main_dfs_dict, peaks_dfs_dict, experiment, organism="fly"
    )


def fly_add_bigWig_data(
    main_dfs_dict: Dict[str, Callable[[], Any]],
    bigWig_data_dict: dict,
    experiment: str,
    res: int,
) -> pd.DataFrame:
    """
    Count statistics (weighted mean, arithmetic mean, minimum and maximum)
    of the bigWig data in both regions of each chromatin loop,
    for each dataframe from the main_dfs_dict dictionary.
    Args:
        main_dfs_dict: dictionary with cell types as keys and load functions of pandas DataFrames with chromatin loops as values.
        bigWig_data_dict: dictionary with cell types as keys and pyBigWig objects as values.
        experiment: name of the experiment
        res: resolution
    Returns:
        dictionary:
            keys: cell types
            values: pandas DataFrames with added columns of bigWig data statistics in both regions of each loop
    """
    if len(bigWig_data_dict) > 1:
        assert (
            len(main_dfs_dict) == 1
        ), "The fly pipeline is only suitable for 1 type of cell."
        cell_name = list(main_dfs_dict.keys())[0]
        for name, path in bigWig_data_dict.items():
            new_experiment = f"{name}_{experiment}"
            main_dfs_dict = add_bigWig_data(
                main_dfs_dict,
                {cell_name: path},
                new_experiment,
                res=res,
                organism="fly",
            )
        return main_dfs_dict
    else:
        return add_bigWig_data(
            main_dfs_dict, bigWig_data_dict, experiment, res=res, organism="fly"
        )


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
    return all_anchors2one_df(dfs_dict, organism="fly")


def fly_all_peaks2one_df(peaks_dict: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    """
    Concatenates all peaks DataFrames into one DataFrame.
    Args:
        peaks_dict: dictionary with cell types as keys and pandas DataFrames as values.
    Returns:
        pandas DataFrame with all peaks.
    """
    return all_peaks2one_df(peaks_dict, organism="fly")


def fly_get_overlaps_with_names(
    anchors_df: pd.DataFrame, peaks_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Get overlapping regions between anchors and peaks.
    Args:
        anchors_df: pandas DataFrame with anchors coordinates.
        peaks_df: pandas DataFrame with peaks coordinates.
    Returns:
        pandas DataFrame with overlapping regions coordinates and names.
    """
    return get_overlaps_with_names(anchors_df, peaks_df)


def fly_getfasta_bedfile(
    df: pd.DataFrame, path_simp_genome: str, path_to_save: str
) -> str:
    """
    Cut sequences from chromosomes using BEDTools for coordinates from the pandas DataFrame.
    Args:
        df: pandas DataFrame with coordinates.
        path_simp_genome: path to fasta file with chromosomes.
        path_to_save: path to save fasta file.
    Returns:
        path_to_save: path to fasta file.
    """
    return getfasta_bedfile(df, path_simp_genome, path_to_save, organism="fly")


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
    return count_motifs(main_dfs_dict, motifs_df, organism="fly")


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


def fly_concat_dfs_from_dict(main_dfs_dict: dict) -> pd.DataFrame:
    """
    Concatenates dataframes from dictionary.
    Args:
        main_dfs_dict: dictionary with cell types as keys and load functions of pandas DataFrames with chromatin loops as values.
        cells_to_use: list of cell types to be used. If empty, all cell types from main_dfs_dict will be used.
    Returns:
        pandas DataFrame with concatenated dataframes from dictionary.
    """
    cells_to_use = list(main_dfs_dict.keys())
    return concat_dfs_from_dict(main_dfs_dict, cells_to_use, organism="fly")
