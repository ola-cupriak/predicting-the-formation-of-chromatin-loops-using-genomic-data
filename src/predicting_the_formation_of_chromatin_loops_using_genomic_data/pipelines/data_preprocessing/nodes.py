import copy
import numpy as np
import math
import pyBigWig
from predicting_the_formation_of_chromatin_loops_using_genomic_data.utils import (
    _dict_partitions,
)
from typing import Any, Callable, Dict
import pandas as pd
from tqdm import tqdm
from Bio import SeqIO
import polars as pl
import pyarrow.parquet as pq
import pybedtools
import subprocess
import warnings
import time

warnings.simplefilter(action="ignore")


def _sort_df(
    df: pd.DataFrame, region_name: str, organism: str = "human"
) -> pd.DataFrame:
    """
    Sort the dataframe by chromosome and position.
    Args:
        df: pandas DataFrame.
        region_name: name of the column with regions.
        organism: name of the organism.
    Returns:
        sorted pandas DataFrame.
    """
    if organism == "human":
        df.loc[df["chr"] == "X", "chr"] = "100"
        df.loc[df["chr"] == "Y", "chr"] = "200"
        try:
            df["chr"] = df["chr"].astype(int)
        except:
            pass
        df = df.sort_values(by=["chr", region_name])
        df.loc[df["chr"] == 100, "chr"] = "X"
        df.loc[df["chr"] == 200, "chr"] = "Y"
        df["chr"] = df["chr"].astype("string")
    else:
        df = df.sort_values(by=["chr", region_name])

    return df


def read_hic(
    partitioned_input: Dict[str, Callable[[], Any]],
    cells2names: Dict[str, dict],
    dataset_name: str,
    r: int,
    cells_to_use: list = None,
) -> Dict[str, pd.DataFrame]:
    """
    Load and modify the dataframes with chromatin loops anotations.
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
    cells_to_use = cells_to_use or []
    dfs_dict = _dict_partitions(partitioned_input)
    cells2names_dataset_dict = cells2names[dataset_name]
    assert set(cells_to_use).issubset(
        set(cells2names_dataset_dict.values())
    ), f"Cell types: {set(cells_to_use)-set(cells2names_dataset_dict.values())} are not in the dataset. Please check data_preprocessing.yml file."
    new_dfs_dict = dict()
    for name, df in dfs_dict.items():
        cell_type = cells2names_dataset_dict[name]
        if cells_to_use and cell_type not in cells_to_use:
            continue
        f1 = lambda x: x["chr1"].split("chr")[-1]
        f2 = lambda x: x["chr2"].split("chr")[-1]
        df["chr1"] = df.apply(f1, axis=1)
        df["chr2"] = df.apply(f2, axis=1)
        assert len(df[df["chr1"] != df["chr2"]]) == 0
        df.rename(columns={"chr1": "chr"}, inplace=True)
        for region in ["x", "y"]:
            df[f"{region}"] = (round((df[f"{region}1"] + df[f"{region}2"]) / 2)).astype(
                int
            )
            df[f"{region}_start"] = (df[f"{region}"] - r).astype(int)
            df[f"{region}_end"] = (df[f"{region}"] + r).astype(int)
        df["len_anchors"] = (df["x2"] - df["x1"]) + (
            df["y2"] - df["y1"]
        )  # add length of both anchors
        df["len_loop"] = df["y2"] - df["x1"]  # add length of loop
        df["cell_type"] = cell_type
        df = df[
            df.columns.intersection(
                ["x", "y", "x_start", "x_end", "y_start", "y_end", "chr", "cell_type"]
            )
        ]
        # sort by chr and x
        df = _sort_df(df, "x")
        # set dtypes
        df = df.astype(
            {
                "x": "int32",
                "y": "int32",
                "x_start": "int32",
                "x_end": "int32",
                "y_start": "int32",
                "y_end": "int32",
                "chr": "string",
                "cell_type": "string",
            }
        )
        new_dfs_dict[cell_type] = df

    return new_dfs_dict


def read_peaks(
    partitioned_input: Dict[str, Callable[[], Any]],
    cells2names: Dict[str, dict],
    dataset_name: str,
    cells_to_use: list,
    organism: str = "human",
) -> Dict[str, pd.DataFrame]:
    """
    Load dataframes with DNase-seq/ChIP-seq peaks and modify the chromosome columns.
    Args:
        partitioned_input: dictionary with partition ids as keys and load functions as values.
        cells2names: dictionary, template: {'dataset_name': {'file_name': 'cell_type'}}
        dataset_name: name of the dataset to select from cells2names.
        cells_to_use: list of cell types to use.
        organism: name of the organism.
    Returns:
        dictionary:
            keys: cell types
            values: pandas DataFrames with DNase-seq/ChIP-seq peaks.
    """
    dfs_dict = _dict_partitions(partitioned_input)
    cells2names_dataset_dict = cells2names[dataset_name]
    assert set(cells_to_use).issubset(
        set(cells2names_dataset_dict.values())
    ), f"Cell types: {set(cells_to_use)-set(cells2names_dataset_dict.values())} are not in the dataset. Please check data_preprocessing.yml file."
    keys_dict = {".".join(key.split(".")[:-1]): key for key in cells2names_dataset_dict}
    new_dfs_dict = dict()
    for name, df in dfs_dict.items():
        df = pl.from_pandas(df)
        cell_type = cells2names_dataset_dict[keys_dict[name]]
        if cells_to_use and cell_type not in cells_to_use:
            continue
        df = df.with_columns(
            df.select(["chr"])
            .apply(lambda x: x[0].split("chr")[-1], return_dtype=pl.Utf8)
            .rename({"apply": "chr"})
        )

        df = df.with_columns(pl.lit(cell_type).cast(pl.Utf8).alias("cell_type"))
        df = df.to_pandas()
        df = _sort_df(df, "start", organism=organism)

        new_dfs_dict[cell_type] = df

    return new_dfs_dict


def read_bigWig(
    partitioned_input: Dict[str, Callable[[], Any]],
    cells2names: Dict[str, dict],
    dataset_name: str,
    cells_to_use: list,
) -> Dict[str, pd.DataFrame]:
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
    input_dict = _dict_partitions(partitioned_input)
    cells2names_dataset_dict = cells2names[dataset_name]
    assert set(cells_to_use).issubset(
        set(cells2names_dataset_dict.values())
    ), f"Cell types: {set(cells_to_use)-set(cells2names_dataset_dict.values())} are not in the dataset. Please check data_preprocessing.yml file."
    keys_dict = {".".join(key.split(".")[:-1]): key for key in cells2names_dataset_dict}
    bigWig_data_dict = dict()
    for name, bigWig_path in input_dict.items():
        try:
            cell_type = cells2names_dataset_dict[keys_dict[name]]
        except:
            continue
        if cells_to_use and cell_type not in cells_to_use:
            continue
        bigWig_data_dict[cell_type] = bigWig_path

    return bigWig_data_dict


def _concat_dfs(dfs_dict: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    """
    Concatenate a dataframes from a dictionary to one dataframe.
    Args:
        dictionary with cell types as keys and pandas DataFrames as values.
    Returns:
        concatenated pandas DataFrame.
    """
    df = pd.concat(list(dfs_dict.values()))

    return df


def _get_overlapping_regions(
    df1: pd.DataFrame,
    df2: pd.DataFrame,
    names: list,
    count: bool = False,
    wa: bool = False,
    wb: bool = False,
) -> pd.DataFrame:
    """
    Get overlapping regions from two bed files.
    Args:
        df1: path to bed file 1
        df2: path to bed file 2
        count: for each entry in A, report the number of hits in B
        wa: write the original entry in A for each overlap.
        wb: write the original entry in B for each overlap.
    Returns:
        pandas DataFrame with overlapping regions or number of overlapping regions.
    """
    bed1 = pybedtools.BedTool.from_dataframe(df1)
    bed2 = pybedtools.BedTool.from_dataframe(df2)

    intersection = bed1.intersect(bed2, c=count, wa=wa, wb=wb)
    intersection = pd.read_table(intersection.fn, names=names)

    return intersection


def _create_pairs_each_with_each_single_df(
    df_pos: pd.DataFrame,
    df_open_chromatin: pd.DataFrame,
    cell_type: str,
    r: int,
    neg_pos_ratio: float,
    random_state: int,
) -> pd.DataFrame:
    """
    Create new anchor pairs from open chromatin regions to generate negative examples.
    Args:
        df_pos: pandas DataFrame with positive anchor pairs.
        df_open_chromatin: pandas DataFrame with open chromatin regions.
        cell_type: cell type to use.
        r: radius of the anchor region.
        neg_pos_ratio: ratio of negative to positive examples.
        random_state: random state.
    Returns:
        pandas DataFrame with positive and negative anchor pairs.
    """
    df_open_chromatin["center"] = (
        df_open_chromatin["start"] + df_open_chromatin["end"]
    ) // 2
    df_open_chromatin = df_open_chromatin[["chr", "center"]]
    chromosomes = pd.unique(df_pos["chr"])
    chromosomes = {k: None for k in chromosomes}
    if neg_pos_ratio:
        chromosomes = {
            k: round(len(df_pos[df_pos["chr"] == k]) * neg_pos_ratio)
            for k in chromosomes
        }
    df_neg = pd.DataFrame()

    for ch in chromosomes.keys():
        df_chr = df_open_chromatin[df_open_chromatin["chr"] == ch]
        df_chr = pl.from_pandas(df_chr)
        if len(df_chr) > 1000:
            df_chr = df_chr.sample(
                n=1000, with_replacement=False, seed=random_state, shuffle=True
            )
        df_chr = df_chr.join(df_chr, how="outer", on="chr", suffix="_y")
        # rename columns in polars
        df_chr = df_chr.rename({"center": "x", "center_y": "y"})
        # filter columns in polars
        df_chr = df_chr.filter(pl.col("x") < pl.col("y"))
        df_chr = df_chr.filter(pl.col("x") + r < pl.col("y") - r)
        # sample with polars
        if neg_pos_ratio:
            df_chr = df_chr.sample(
                n=chromosomes[ch],
                with_replacement=False,
                seed=random_state,
                shuffle=True,
            )
        # convert back to pandas
        df_chr = df_chr.to_pandas()
        df_neg = pd.concat([df_neg, df_chr])

    df_neg["cell_type"] = cell_type
    df_neg["label"] = 0
    df_neg["x_start"] = df_neg["x"] - r
    df_neg["x_end"] = df_neg["x"] + r
    df_neg["y_start"] = df_neg["y"] - r
    df_neg["y_end"] = df_neg["y"] + r

    df_pairs = pd.concat([df_pos, df_neg])
    df_pairs.reset_index(drop=True, inplace=True)

    df_pairs = df_pairs.astype(
        {
            "x": "int32",
            "y": "int32",
            "x_start": "int32",
            "x_end": "int32",
            "y_start": "int32",
            "y_end": "int32",
            "label": "int16",
        }
    )

    return df_pairs


def _x_and_y_anchors2one_col(
    df: pd.DataFrame,
    x_cols2use: list = None,
    y_cols2use: list = None,
    new_cols: list = None,
    sortby: str = "start",
    organism: str = "human",
) -> pd.DataFrame:
    """
    Combines the columns describing x regions and the columns describing y regions
    into one set of columns describing all regions.
    (columns: x_chr, x_start, x_end, y_chr, y_start, y_end -> columns: chr, start, end)
    Args:
        df: DataFrame with columns describing x regions and columns describing y regions.
        x_cols2use: columns describing x regions.
        y_cols2use: columns describing y regions.
        new_cols: new column names.
        sortby: sort by this column.
        organism: name of the organism.
    Returns:
        pandas DataFrame with one set of columns describing all regions.
    """
    x_cols2use = x_cols2use or ["chr", "x_start", "x_end"]
    y_cols2use = y_cols2use or ["chr", "y_start", "y_end"]
    new_cols = new_cols or ["chr", "start", "end"]

    df_part1 = df[x_cols2use]
    df_part2 = df[y_cols2use]
    # rename columns
    df_part1.columns = new_cols
    df_part2.columns = new_cols
    # combine two columns into one
    anchors_df = pd.concat([df_part1, df_part2], axis=0)
    # remove duplicate rows
    anchors_df = anchors_df.drop_duplicates()
    # sort by chr and region
    anchors_df = _sort_df(anchors_df, sortby, organism=organism)
    # reset index
    anchors_df = anchors_df.reset_index(drop=True)
    anchors_df = anchors_df[new_cols]

    return anchors_df


def _create_new_anchors_pairs(
    df: pd.DataFrame, organism: str = "human"
) -> pd.DataFrame:
    """
    Join all values from one column (x) with all values from second column (y).
    Args:
        df: pandas DataFrame with columns describing x regions and columns describing y regions.
        organism: name of the organism.
    Returns:
        pandas DataFrame with new pairs of anchors.
    """
    df_new_pairs = _x_and_y_anchors2one_col(
        df, ["chr", "x"], ["chr", "y"], ["chr", "anchor"], "anchor", organism=organism
    )
    df_new_pairs = pl.from_pandas(df_new_pairs)
    df_new_pairs = df_new_pairs.join(
        df_new_pairs, how="outer", on="chr", suffix="_y"
    ).filter(pl.col("anchor") < pl.col("anchor_y"))
    df_new_pairs = df_new_pairs.rename({"anchor": "x", "anchor_y": "y"})
    df_new_pairs = df_new_pairs.to_pandas()

    return df_new_pairs


def _remove_duplicated_pairs(
    df_pos: pd.DataFrame, df_neg: pd.DataFrame
) -> pd.DataFrame:
    """
    Remove pairs from neg_df that are present in pos_df
    Args:
        df_pos: pandas DataFrame with positive examples.
        df_neg: pandas DataFrame with negative examples.
    Returns:
        pandas DataFrame with negative examples that are not present in positive examples.
    """
    df_pos = df_pos[["chr", "x", "y"]]
    df_neg = df_neg[["chr", "x", "y"]]
    df_neg = df_neg.merge(df_pos, how="left", on=["chr", "x", "y"], indicator=True)
    df_neg = df_neg[df_neg["_merge"] == "left_only"]
    df_neg.drop(columns="_merge", inplace=True)

    return df_neg


def _get_negatives_by_new_anchors_pairing(
    df: pd.DataFrame,
    cell_type: str,
    r: int,
    neg_pos_ratio: float,
    random_state: int,
    df_len: int,
    organism: str = "human",
    df_neg: pd.DataFrame = None,
) -> pd.DataFrame:
    """
    Create negative examples by pairing anchors derived from the same chromosome.
    Args:
        df: pandas DataFrame with positive examples.
        cell_type: cell type of the positive examples.
        r: radius of the region around the anchor.
        neg_pos_ratio: ratio of negative to positive examples.
        random_state: random state.
        df_len: length of the positive examples DataFrame.
        organism: name of the organism.
    Returns:
        pandas DataFrame with positive and negative examples.
    """
    if df_neg is None:
        df_neg = _create_new_anchors_pairs(df, organism=organism)
        df_neg = _remove_duplicated_pairs(df, df_neg)
    df_neg["x_start"] = df_neg["x"] - r
    df_neg["x_end"] = df_neg["x"] + r
    df_neg["y_start"] = df_neg["y"] - r
    df_neg["y_end"] = df_neg["y"] + r
    df_neg = df_neg[df_neg["x_end"] < df_neg["y_start"]]
    df_neg["cell_type"] = cell_type
    df_neg["label"] = 0
    df_neg = df_neg[
        ["chr", "x", "x_start", "x_end", "y", "y_start", "y_end", "cell_type", "label"]
    ]
    # randomy sample negative examples WITHOUT REPEATS to get the desired neg_pos_retio
    if neg_pos_ratio:
        df_neg = df_neg.sample(
            n=round(int(neg_pos_ratio * df_len)),
            random_state=random_state,
            replace=False,
        )
    # merge positive and negative examples
    df = pd.concat([df, df_neg], axis=0)
    # sort by chr and region
    df = _sort_df(df, "x_start", organism=organism)
    df = df.reset_index(drop=True)
    df = df.astype(
        {
            "x": "int32",
            "y": "int32",
            "x_start": "int32",
            "x_end": "int32",
            "y_start": "int32",
            "y_end": "int32",
            "label": "int16",
        }
    )

    return df


def _add_distances(df: pd.DataFrame) -> pd.DataFrame:
    """
    Add column with distances between x and y anchors.
    Args:
        df: pandas DataFrame with anchors.
    Returns:
        pandas DataFrame with distances between x and y anchors.
    """
    df["distance"] = df["y"] - df["x"]
    df = df.astype({"distance": "int32"})

    return df


def _remove_negative_anchors_overlapping_positives(df: pd.DataFrame):
    """
    Remove negative anchors that overlap with positive anchors.
    Args:
        df: pandas DataFrame.
    Returns:
        pandas DataFrame with removed negative anchors that overlap with positive anchors.
    """
    # create new df with positive and negative anchors (both x and y) in the same columns
    df_pos = _x_and_y_anchors2one_col(df[df["label"] == 1])
    df_pos["id"] = [i for i in range(len(df_pos))]
    df_neg = _x_and_y_anchors2one_col(df[df["label"] == 0])
    df_neg["id"] = [i for i in range(len(df_neg))]
    # find negatives overlapping with positives
    df_overlap = _get_overlapping_regions(
        df_pos,
        df_neg,
        names=["chr1", "start1", "end1", "id_pos", "chr2", "start2", "end2", "id_neg"],
        wa=True,
        wb=True,
    )
    # drop negatives overlapping with positives from df_neg
    df_neg = df_neg.loc[~df_neg["id"].isin(df_overlap["id_neg"]), :]
    df_neg.drop(["id"], axis=1, inplace=True)
    # add x and y anchors to df_neg from df
    df_neg["center"] = (round((df_neg[f"start"] + df_neg[f"end"]) / 2)).astype(int)
    return df_neg


def _get_negatives_by_random_pairing_other_cell_types(
    df: pd.DataFrame, neg_pos_ratio: int, random_state: int, cell_type: str, r: int
):
    df_neg = _remove_negative_anchors_overlapping_positives(df)
    # random choose half of the negative examples
    df_neg = df_neg.sample(
        frac=0.5,
        random_state=random_state,
        replace=False,
    )
    df = df[df["label"] == 1]
    # Create pairs of anchors
    df_neg = df_neg[["chr", "center"]]
    df_neg = pl.from_pandas(df_neg)
    df_neg = df_neg.join(df_neg, how="outer", on="chr", suffix="_y").filter(
        pl.col("center") < pl.col("center_y")
    )
    df_neg = df_neg.rename({"center": "x", "center_y": "y"})
    df_neg = df_neg.to_pandas()
    df_neg = df_neg[df_neg["x"] < df_neg["y"]]
    # Randomly sample negative examples WITHOUT REPEATS to get the desired neg_pos_retio
    pos_num = len(df)
    df = _get_negatives_by_new_anchors_pairing(
        df, cell_type, r, neg_pos_ratio, random_state, pos_num, df_neg=df_neg
    )

    return df


def add_labels(
    dfs_dict: Dict[str, pd.DataFrame],
    negtype: str,
    peaks_dict: Dict[str, pd.DataFrame],
    r: int,
    neg_pos_ratio: float,
    random_state: int,
    organism: str = "human",
) -> Dict[str, pd.DataFrame]:
    """
    Add labels to the dataframes depending on the cell type and type of model to train.
    Within model: 1 if cell type is the same as the cell type of the loop, 0 otherwise.
    Across model: negative sampling involving the use of open regions of chromatin that
                do not overlap with any positive example
    Args:
        dfs_dict: dictionary with cell types as keys and pandas DataFrames as values.
        negtype: type of negative sampling.
        peaks_dict: dictionary with cell types as keys and pandas DataFrames with peaks as values.
        r: radius of the region around the anchor.
        neg_pos_ratio: ratio of negative to positive examples.
        random_state: random state.
        organism: name of the organism.
    Returns:
        dictionary:
            keys: cell types
            values: pandas DataFrames with labels added.
    """
    assert negtype in [
        "anchors_from_other_cell_types",
        "new_pairs_of_anchors",
        "open_chromatin_regions",
        "new_pairs_and_open_chromatin",
    ], f"Negative sampling type should be either 'anchors_from_other_cell_types', 'new_pairs_of_anchors' or 'open_chromatin_regions, but got {type}."
    if negtype == "anchors_from_other_cell_types":
        df = _concat_dfs(dfs_dict)
        df = _sort_df(df, "x", organism=organism)

    for name, cell_df in dfs_dict.items():
        print("Creating negative examples for cell type", name, "...")

        if negtype == "anchors_from_other_cell_types":
            f = lambda x: 1 if x["cell_type"] == name else 0
            df["label"] = df.apply(f, axis=1)
            df["label"] = df["label"].astype("int16")
            dfs_dict[name] = _get_negatives_by_random_pairing_other_cell_types(
                copy.deepcopy(df), neg_pos_ratio, random_state, name, r
            )
        elif negtype == "new_pairs_of_anchors":
            cell_df["label"] = 1
            df_len = len(cell_df)
            dfs_dict[name] = _get_negatives_by_new_anchors_pairing(
                cell_df, name, r, neg_pos_ratio, random_state, df_len, organism=organism
            )
        elif negtype == "open_chromatin_regions":
            assert (
                peaks_dict != None
            ), "peaks_dict should be provided for this type of negative sampling"
            cell_df["label"] = 1
            dfs_dict[name] = _create_pairs_each_with_each_single_df(
                cell_df, peaks_dict[name], name, r, neg_pos_ratio, random_state
            )
        else:
            assert (
                peaks_dict != None
            ), "peaks_dict should be provided for this type of negative sampling"
            cell_df["label"] = 1
            df_len = len(cell_df)
            with_neg = _create_pairs_each_with_each_single_df(
                cell_df, peaks_dict[name], name, r, neg_pos_ratio / 2, random_state
            )
            with_neg = _get_negatives_by_new_anchors_pairing(
                with_neg,
                name,
                r,
                neg_pos_ratio / 2,
                random_state,
                df_len,
                organism=organism,
            )
            dfs_dict[name] = with_neg

    return dfs_dict


def _find_the_closest_peak(
    main_df: pd.DataFrame, peak_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Find the distance of the closest peak from the centre of each region.
    Args:
        main_df: pandas DataFrame with regions.
        peak_df: pandas DataFrame with peaks.
    Returns:
        pandas DataFrame with the distances of the closest peaks to the centre of each region.
    """
    # Remove peaks on chromosomes that are not in the main_df
    chromosomes = main_df["chr"].unique()
    peak_df = peak_df.loc[peak_df["chr"].isin(chromosomes), :]

    assert (
        len(main_df.loc[main_df.iloc[:, 1] != main_df.iloc[:, 2], :]) == 0
    ), "The region columns are not equal"
    bed1 = pybedtools.BedTool.from_dataframe(main_df)
    bed2 = pybedtools.BedTool.from_dataframe(peak_df)

    closest = bed1.closest(bed2, d=True, t="first")
    closest = pd.read_table(
        closest.fn,
        names=[
            "chr",
            "centre",
            "centre_dup",
            "second_reg",
            "cell_type",
            "chr_peak",
            "start_peak",
            "end_peak",
            "cell_type_peak",
            "distance",
        ],
    )
    closest = closest.loc[:, ["chr", "centre", "second_reg", "distance", "cell_type"]]

    return closest


def _count_peaks_single_df(
    main_df: pd.DataFrame, peaks_df: pd.DataFrame, experiment: str
) -> pd.DataFrame:
    """
    Count the number of experiment peaks in both regions of each chromatin loop.
    Args:
        main_df: pandas DataFrame with chromatin loops.
        peaks_df: pandas DataFrame with experiment peaks.
        experiment: name of the experiment.
    Returns:
        pandas DataFrame with added columns of numbers of experiment peaks in both regions of each loop
    """
    for region in ["x", "y"]:
        second_reg = "y" if region == "x" else "x"
        peak_counts = main_df.loc[
            :, ["chr", f"{region}_start", f"{region}_end", f"{second_reg}", "cell_type"]
        ]
        peak_counts = _get_overlapping_regions(
            peak_counts,
            peaks_df,
            names=["chr", "start", "end", "second_reg", "cell_type", "count"],
            count=True,
        )
        peak_counts["chr"] = peak_counts["chr"].astype("str")
        peak_counts.rename(
            columns={
                "count": f"{region}_{experiment}_counts",
                "start": f"{region}_start",
                "end": f"{region}_end",
                "second_reg": f"{second_reg}",
            },
            inplace=True,
        )

        assert len(peak_counts) == len(
            main_df
        ), "Length of the main_df and peak_counts are not equal"

        main_len_before = len(main_df)
        main_df = main_df.merge(peak_counts)

        assert (
            len(main_df) == main_len_before
        ), "Length of the main_df changed after merging"

    to_change_dtype = [f"x_{experiment}_counts", f"y_{experiment}_counts"]
    main_df.loc[:, to_change_dtype] = main_df.loc[:, to_change_dtype].astype("int16")

    return main_df


def _find_the_closest_peaks_single_df(
    main_df: pd.DataFrame,
    peaks_df: pd.DataFrame,
    experiment: str,
    organism: str = "human",
) -> pd.DataFrame:
    """
    Find the distance of the closest peak of the experiment from the center of each anchor from the chromatin loop.
    Args:
        main_df: pandas DataFrame with chromatin loops.
        peaks_df: pandas DataFrame with experiment peaks.
        experiment: name of the experiment.
        organism: name of the organism.
    Returns:
        pandas DataFrame with added columns of distances of the closest peak of the experiment
        from the center of each anchor from the chromatin loop.
    """
    for region in ["y", "x"]:
        second_reg = "y" if region == "x" else "x"
        distances = main_df.loc[:, ["chr", f"{region}", f"{second_reg}", "cell_type"]]
        distances = _sort_df(distances, region, organism=organism)
        distances = distances.loc[
            :, ["chr", f"{region}", f"{region}", f"{second_reg}", "cell_type"]
        ]
        distances = _find_the_closest_peak(distances, peaks_df)
        distances["chr"] = distances["chr"].astype("str")
        distances.rename(
            columns={
                "distance": f"{region}_{experiment}_distance",
                "centre": f"{region}",
                "second_reg": f"{second_reg}",
            },
            inplace=True,
        )

        assert len(distances) == len(
            main_df
        ), "Length of the main_df and distances are not equal"

        main_len_before = len(main_df)
        main_df = main_df.merge(distances)

        assert (
            len(main_df) == main_len_before
        ), "Length of the main_df changed after merging"

    to_change_dtype = [f"x_{experiment}_distance", f"y_{experiment}_distance"]
    main_df.loc[:, to_change_dtype] = main_df.loc[:, to_change_dtype].astype("int32")

    return main_df


def count_peaks_and_distances(
    main_dfs_dict: Dict[str, Callable[[], Any]],
    peaks_dfs_dict: Dict[str, pd.DataFrame],
    experiment: str,
    organism: str = "human",
) -> Dict[str, pd.DataFrame]:
    """
    Add a distance beetwen the anchors,
    count the number of peaks in both regions of each chromatin loop
    and find the distance of the closest peak from the center of each anchor from the chromatin loop,
    for each dataframe from the main_dfs_dict dictionary.
    Args:
        main_dfs_dict: dictionary with cell types as keys and load functions of pandas DataFrames with chromatin loops as values.
        peaks_dfs_dict: dictionary with cell types as keys and pandas DataFrames with experiment peaks as values.
        experiment: name of the experiment.
        organism: name of the organism.
    Returns:
        dictionary:
            keys: cell types
            values: pandas DataFrame with added columns of numbers of experiment peaks in both regions of each chromatin loop
                    and columns of distances of the closest peak of the experiment from the center of each anchor from the chromatin loop.
    """
    print(f"Adding peaks counts and distances from anchor centers for {experiment}...")
    # if main_dfs_dict values are not DataFrames, load them
    if not isinstance(list(main_dfs_dict.values())[0], pd.DataFrame):
        main_dfs_dict = _dict_partitions(main_dfs_dict)
    # add distances
    for cell_type, df in main_dfs_dict.items():
        main_dfs_dict[cell_type] = _add_distances(df)
    for peaks_name, peaks_df in peaks_dfs_dict.items():
        print(f"...for {peaks_name} cell...")
        for main_name, main_df in main_dfs_dict.items():
            if main_name == peaks_name:
                main_df = _count_peaks_single_df(main_df, peaks_df, experiment)
                main_df = _find_the_closest_peaks_single_df(
                    main_df, peaks_df, experiment, organism=organism
                )
                assert (
                    len(main_df[main_df[f"x_{experiment}_distance"] < 0]) == 0
                ), f"Negative distances found in column x_{experiment}_distance for cell {main_name}"
                assert (
                    len(main_df[main_df[f"y_{experiment}_distance"] < 0]) == 0
                ), f"Negative distances found in column y_{experiment}_distance for cell {main_name}"
                main_dfs_dict[main_name] = main_df
    print("Done!")

    return main_dfs_dict


def _get_triangle_kernel(kerlen: int) -> list:
    """
    Get a triangle kernel of length kerlen.
    Args:
        kerlen: length of the kernel.
    Returns:
        list with the triangle kernel.
    """
    r = np.arange(kerlen)
    kernel1d = (kerlen + 1 - np.abs(r - r[::-1])) / 2
    kernel1d /= kernel1d.sum()

    return kernel1d


def _calculate_weighted_mean(distribution: list):
    """
    Calculate weighted mean of the distribution.
    Args:
        distribution: list with values of the distribution.
    Returns:
        weighted mean of the distribution.
    """
    weights = _get_triangle_kernel(len(distribution))
    return sum([distribution[i] * weights[i] for i in range(len(distribution))]) / sum(
        weights
    )


def _replace_nans_with_zeros(distribution: list) -> list:
    """
    Replace NaNs with zeros in the distribution.
    Args:
        distribution: list with values of the distribution.
    Returns:
        list with zeros instead of NaNs.
    """
    return [0 if np.isnan(x) else x for x in distribution]


def _get_stats_single_row(bigwig_obj, chromosome, start, end, res) -> tuple:
    """
    Get statistics of the distribution.
    Args:
        distribution: list with values of the distribution.
    Returns:
        tuple with mean, weighted mean, minimum and maximum of the distribution.
    """
    try:
        distribution = bigwig_obj.values(chromosome, start, end)
    except:
        chromosome = chromosome.split("chr")[-1]
        distribution = bigwig_obj.values(chromosome, start, end)
    distribution = _replace_nans_with_zeros(distribution)
    mean = sum(distribution) / len(distribution)
    weighted_mean = _calculate_weighted_mean(distribution)
    # add resolution to minimum and maximum
    assert (
        len(distribution) % res == 0
    ), "Length of the region is not divisible by resolution"
    n = len(distribution) // res
    means = []
    for i in range(n):
        means.append(sum(distribution[i * res : (i + 1) * res]) / res)
    maximum = max(means)
    minimum = min(means)

    return mean, weighted_mean, minimum, maximum


def _get_stats_all_anchors(main_df: pl.DataFrame, bigWig_obj, res: int) -> pd.DataFrame:
    """
    Get statistics of the distribution for all anchors.
    Args:
        main_df: polars DataFrame with anchors.
        bigWig_obj: pyBigWig object with bigwig file.
        res: resolution
    Returns:
        polars DataFrame with mean, weighted mean, minimum and maximum of the distribution for all anchors.
    """
    df_x = main_df.select(["chr", "x_start", "x_end"])
    df_y = main_df.select(["chr", "y_start", "y_end"])
    df_y = df_y.rename({"y_start": "x_start", "y_end": "x_end"})
    # polars concat
    df = pl.concat([df_x, df_y])
    df = df.unique(subset=["chr", "x_start", "x_end"])
    df = df.with_columns(
        df.select(["chr", "x_start", "x_end"])
        .apply(lambda x: _get_stats_single_row(bigWig_obj, x[0], x[1], x[2], res))
        .rename(
            {
                "column_0": "mean",
                "column_1": "weighted_mean",
                "column_2": "min",
                "column_3": "max",
            }
        )
    )

    return df


def _regions_mean_of_means(main_df: pd.DataFrame, experiment: str) -> float:
    """
    Calculate mean of means of the experiment in both regions of each chromatin loop.
    Args:
        main_df: pandas DataFrame with chromatin loops.
        experiment: name of the experiment.
    Returns:
        mean of means of the experiment in both regions of each chromatin loop.
    """
    df_x = main_df.loc[:, ["chr", "x", f"x_{experiment}_mean"]]
    df_y = main_df.loc[:, ["chr", "y", f"y_{experiment}_mean"]]
    df_y.rename(
        columns={"y": "x", f"y_{experiment}_mean": f"x_{experiment}_mean"}, inplace=True
    )
    df = pd.concat([df_x, df_y], axis=0)
    df = df.drop_duplicates(subset=["chr", "x"])
    mean = df[f"x_{experiment}_mean"].mean()

    return mean


def _add_bigWig_data_single_df(
    main_df: pd.DataFrame, bigWig_path, experiment: str, organism: str, res: int
) -> pd.DataFrame:
    """
    Count statistics (weighted mean, arithmetic mean, minimum and maximum)
    of the bigWig data in both regions of each chromatin loop.
    Args:
        main_df: pandas DataFrame with chromatin loops.
        bigWig_obj: pyBigWig object with experiment peaks
        experiment: name of the experiment.
        organism: name of the organism.
        res: resolution
    Returns:
        pandas DataFrame with added columns of bigWig data statistics in both regions of each loop
    """
    bigWig_obj = pyBigWig.open(bigWig_path)
    regions = ["x", "y"]

    main_df["chr"] = main_df["chr"].apply(lambda x: "chr" + x)

    main_df = pl.from_pandas(main_df)
    df_with_stats = _get_stats_all_anchors(main_df, bigWig_obj, res)
    for region in regions:
        if region == "y":
            df_with_stats = df_with_stats.rename(
                {"x_start": "y_start", "x_end": "y_end"}
            )

        main_df = main_df.join(
            df_with_stats, on=["chr", f"{region}_start", f"{region}_end"], how="left"
        )
        main_df = main_df.rename(
            {
                "mean": f"{region}_{experiment}_mean",
                "weighted_mean": f"{region}_{experiment}_weighted_mean",
                "min": f"{region}_{experiment}_min",
                "max": f"{region}_{experiment}_max",
            }
        )

        main_df = main_df.with_columns(
            pl.col(f"{region}_{experiment}_mean").cast(pl.Float32),
            pl.col(f"{region}_{experiment}_weighted_mean").cast(pl.Float32),
            pl.col(f"{region}_{experiment}_min").cast(pl.Float32),
            pl.col(f"{region}_{experiment}_max").cast(pl.Float32),
        )

    main_df = main_df.to_pandas()
    mean = _regions_mean_of_means(main_df, experiment)
    # divide all created columns by mean of means
    for region in regions:
        main_df[f"{region}_{experiment}_mean"] = (
            main_df[f"{region}_{experiment}_mean"] / mean
        )
        main_df[f"{region}_{experiment}_weighted_mean"] = (
            main_df[f"{region}_{experiment}_weighted_mean"] / mean
        )
        main_df[f"{region}_{experiment}_min"] = (
            main_df[f"{region}_{experiment}_min"] / mean
        )
        main_df[f"{region}_{experiment}_max"] = (
            main_df[f"{region}_{experiment}_max"] / mean
        )

    main_df["chr"] = main_df["chr"].apply(lambda x: x.replace("chr", ""))

    assert main_df.isna().any().sum() == 0, "There are NaN values in the dataframe"

    return main_df


def add_bigWig_data(
    main_dfs_dict: Dict[str, Callable[[], Any]],
    bigWig_data_dict: dict,
    experiment: str,
    res: int = 20,
    organism: str = "human",
) -> Dict[str, pd.DataFrame]:
    """
    Count statistics (weighted mean, arithmetic mean, minimum and maximum)
    of the bigWig data in both regions of each chromatin loop,
    for each dataframe from the main_dfs_dict dictionary.
    Args:
        main_dfs_dict: dictionary with cell types as keys and load functions of pandas DataFrames with chromatin loops as values.
        bigWig_data_dict: dictionary with cell types as keys and pyBigWig objects as values.
        experiment: name of the experiment.
        organism: name of the organism.
        res: resolution
    Returns:
        dictionary:
            keys: cell types
            values: pandas DataFrames with added columns of bigWig data statistics in both regions of each loop
    """
    print(f"Adding {experiment} data...")
    # if main_dfs_dict values are not DataFrames, load them
    if not isinstance(list(main_dfs_dict.values())[0], pd.DataFrame):
        main_dfs_dict = _dict_partitions(main_dfs_dict)
    for bigWig_name, bigWig_path in bigWig_data_dict.items():
        print(f"...for {bigWig_name} cell...")
        for main_name, main_df in main_dfs_dict.items():
            if main_name == bigWig_name:
                main_df = _add_bigWig_data_single_df(
                    main_df, bigWig_path, experiment, organism, res
                )
                main_dfs_dict[main_name] = main_df
    print("Done!")

    return main_dfs_dict


def all_anchors2one_df(
    dfs_dict: Dict[str, pd.DataFrame], organism: str = "human"
) -> pd.DataFrame:
    """
    Concatenates all anchors DataFrames from dfs_dict into one DataFrame and
    combines the columns describing x regions and the columns describing y regions
    into one set of columns describing all regions.
    (columns: x_chr, x_start, x_end, y_chr, y_start, y_end -> columns: chr, start, end)
    Args:
        dfs_dict: dictionary with cell types as keys and pandas DataFrames as values.
        organism: name of the organism.
    Returns:
        pandas DataFrame with one set of columns describing all regions.
    """
    # if dfs_dict values are not DataFrames, load them
    if not isinstance(list(dfs_dict.values())[0], pd.DataFrame):
        dfs_dict = _dict_partitions(dfs_dict)

    df_with_2_regions = _concat_dfs(dfs_dict)
    anchors_df = _x_and_y_anchors2one_col(df_with_2_regions, organism=organism)

    return anchors_df


def all_peaks2one_df(
    peaks_dict: Dict[str, pd.DataFrame], organism: str = "human"
) -> pd.DataFrame:
    """
    Concatenates all peaks DataFrames into one DataFrame.
    Args:
        peaks_dict: dictionary with cell types as keys and pandas DataFrames as values.
        organism: name of the organism.
    Returns:
        pandas DataFrame with all peaks.
    """
    df = _concat_dfs(peaks_dict)
    # sort by chr and region
    df = _sort_df(df, "start", organism=organism)
    df = df.reset_index(drop=True)

    return df


def get_overlaps_with_names(
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
    intersection_wa_wb = _get_overlapping_regions(
        anchors_df,
        peaks_df,
        names=[
            "anchor_chr",
            "anchor_start",
            "anchor_end",
            "peak_chr",
            "peak_start",
            "peak_end",
            "cell_type",
        ],
        wa=True,
        wb=True,
    )
    intersection = _get_overlapping_regions(
        anchors_df, peaks_df, names=["chr", "start", "end"]
    )

    joined_intersection = pd.merge(
        intersection_wa_wb, intersection, left_index=True, right_index=True
    )
    joined_intersection["name"] = joined_intersection.apply(
        lambda x: f"{x['chr']}:{x['anchor_start']}:{x['anchor_end']}:{x['cell_type']}",
        axis=1,
    )

    return joined_intersection[["chr", "start", "end", "name"]]


def getfasta_bedfile(
    df: pd.DataFrame, path_simp_genome: str, path_to_save: str, organism: str = "human"
) -> str:
    """
    Cut sequences from chromosomes using BEDTools for coordinates from the pandas DataFrame.
    Args:
        df: pandas DataFrame with coordinates.
        path_simp_genome: path to fasta file with chromosomes.
        path_to_save: path to save fasta file.
        organism: name of the organism.
    Returns:
        path_to_save: path to fasta file.
    """
    fasta = pybedtools.BedTool(path_simp_genome)
    df_to_search = df
    df_to_search["chr"] = "chr" + df_to_search["chr"].astype(str)
    bed = pybedtools.BedTool.from_dataframe(df_to_search)
    fasta_bed = bed.sequence(fi=fasta, nameOnly=True)

    records = []
    for record in SeqIO.parse(fasta_bed.seqfn, "fasta"):
        records.append(record)

    SeqIO.write(records, path_to_save, "fasta")

    return path_to_save


def _modify_fimo_output(line: str, i: int):
    """
    Modifies each line of output of the FIMO tool.
    Saves motif IDs and corresponding names to a file.
    Args:
        lines: list of lines from the output file.
    Returns:
        list of modified lines.
    """
    # Change header to tsv format
    if i == 0:
        line = line.replace("sequence name", "chr\tanchor_start\tanchor_end\tcell_type")
        line = line.replace("#pattern name", "motif_id")
    else:
        # Get the information you need
        elements = line.split("\t")
        motif_id, name, _, _, strand, _, _, _, _ = elements
        # Save motif_id and motif_alt_id to the dictionary
        # Change name format to tsv format
        name = name.replace(":", "\t")
        # Add strand information to the motif_id
        if strand == "+":
            motif_id += "_f"
        else:
            motif_id += "_r"
        # Update the list of lines
        elements[0] = motif_id
        elements[1] = name
        line = "\t".join(elements)

    return line


def find_motifs(path_motifs: str, path_fasta: list) -> pd.DataFrame:
    """
    Finds motif occurances in fasta sequences.
    Args:
        path_motifs: path to file with motifs in jaspar format.
        path_fasta: path to fasta file with sequences.
    Returns:
        pandas DataFrame with counts of each motif found in fasta sequences.
    """
    # Change jaspar format to meme format
    path_for_meme = path_motifs.replace(".txt", ".meme")
    subprocess.run(f"jaspar2meme -bundle {path_motifs} > {path_for_meme}", shell=True)
    # Find motif occurances in fasta sequences
    start = time.time()
    print("Finding motifs...")
    subprocess.run(
        f"fimo --text {path_for_meme} {path_fasta} > data/temp/temp.csv", shell=True
    )
    print(f"Finding motifs took {time.time() - start} seconds")
    # Read output
    lines = []
    with open("data/temp/temp.csv", "r") as f:
        with open("data/temp/temp2.csv", "w") as f2:
            pass
        for i, line in tqdm(enumerate(f)):
            lines.append(_modify_fimo_output(line, i))
            if i % 1000000 == 0:
                print(i)
                # Save modified output to temporary file
                with open("data/temp/temp2.csv", "a") as f2:
                    f2.writelines(lines)
                    del lines
                    lines = []
        with open("data/temp/temp2.csv", "a") as f2:
            f2.writelines(lines)
            del lines
    # Read temporary file as pandas DataFrame
    dtypes = {
        "chr": "string",
        "anchor_start": "int32",
        "anchor_end": "int32",
        "motif_id": "string",
        "cell_type": "string",
        "score": "float32",
    }
    df = pd.read_csv(
        "data/temp/temp2.csv", sep="\t", dtype=dtypes, usecols=list(dtypes.keys())
    )
    df.rename(columns={"anchor_start": "start", "anchor_end": "end"}, inplace=True)
    subprocess.run("rm data/temp/temp.csv", shell=True)
    # subprocess.run('rm data/temp/temp2.csv', shell=True)
    df = df[["chr", "start", "end", "motif_id", "cell_type", "score"]]
    len_before = len(df)
    # Count motif occurences
    df = pd.DataFrame(
        df.groupby(
            ["chr", "start", "end", "cell_type", "motif_id"], observed=True
        ).size(),
        columns=["count"],
    )
    df["count"] = df["count"].astype("int16")
    df = df.unstack("motif_id", fill_value=0)
    df.columns = ["_".join(x) for x in df.columns if x[0] == "count"]
    df.reset_index(inplace=True)
    df.rename(
        columns={name: name.replace("count_", "") for name in df.columns}, inplace=True
    )

    assert (
        sum(df.iloc[:, 4:].sum(axis=1)) == len_before
    ), "Something went wrong with counting motifs"

    return df


def _reverse_names(colnames: list):
    """
    Create dictionary with swapped suffixes.
    Columns with suffix '_f' will be swapped to columns with suffix '_r' and vice versa.
    Args:
        colnames: list of column names.
    Returns:
        dictionary with swapped suffixes.
    """
    to_reverse_f = [name for name in colnames if "_f" in name]
    to_reverse_r = [name for name in colnames if "_r" in name]
    to_reverse_f = {name: name.replace("_f", "_r") for name in to_reverse_f}
    to_reverse_r = {name: name.replace("_r", "_f") for name in to_reverse_r}

    to_reverse = {**to_reverse_f, **to_reverse_r}

    return to_reverse


def _count_motifs_single_df(
    main_df: pd.DataFrame,
    motifs_df: pd.DataFrame,
    cell_type: str,
    organism: str = "human",
) -> pd.DataFrame:
    """
    Counts motif occurances in both regions of each chromatin loop.
    Args:
        main_df: pandas DataFrame with chromatin loops.
        motifs_df: pandas DataFrame with motif occurrences found.
        cell_type: name of the cell type.
        organism: name of the organism.
    Returns:
        pandas DataFrame with added columns of each motif counts in both regions of each chromatin loop.
    """
    # change dtypes in main_df
    main_df[["chr", "cell_type"]] = main_df[["chr", "cell_type"]].astype("string")

    # add id column
    main_df["id"] = [i for i in range(len(main_df))]
    main_df_new = copy.deepcopy(main_df)
    columns_no = len(main_df_new.columns)

    # reverse _f/_r suffix in names of columns with motifs
    motifs_colnames = list(motifs_df.columns)
    motifs_colnames.remove("chr")
    to_reverse = _reverse_names(motifs_colnames)

    for region in ["x", "y"]:
        if region == "y":
            motifs_df = motifs_df.rename(columns=to_reverse)
        motifs_df_renamed = motifs_df.rename(
            columns={name: f"{region}_{name}" for name in motifs_colnames}
        )
        main_df_new = main_df_new.merge(
            motifs_df_renamed,
            on=["chr", f"{region}_start", f"{region}_end"],
            how="inner",
        )
    # add rows with missing ids
    motifs_list = list(main_df_new.columns)[columns_no:]
    to_add = main_df[~main_df["id"].isin(main_df_new["id"])]
    if len(to_add):
        for m in motifs_list:
            to_add.loc[:, m] = pd.Series([0] * len(main_df), dtype="int16")
        main_df_new = pd.concat([main_df_new, to_add], ignore_index=True)
    main_df_new = _sort_df(main_df_new, "x_start", organism=organism)
    main_df_new = main_df_new.loc[:, main_df_new.columns != "id"]

    # dtypes = {'chr': "string", 'anchor_start': "int32", 'anchor_end': "int32", 'motif_id': "string", 'cell_type': "string", "score": "float32"}
    # df_found_motifs = pd.read_csv('data/temp/temp2.csv', sep='\t', dtype=dtypes, usecols=list(dtypes.keys()))
    # df_found_motifs = df_found_motifs[df_found_motifs['cell_type'] == cell_type]
    # df_found_motifs.rename(columns={'anchor_start': 'start', 'anchor_end': 'end'}, inplace=True)
    # df_found_motifs['motif_id_without_orientation'] = [motif_id.split('_')[0] for motif_id in df_found_motifs['motif_id']]
    # idx = df_found_motifs.groupby(['motif_id_without_orientation', 'chr', 'start', 'end'])['score'].transform(max) == df_found_motifs['score']
    # df_found_motifs = df_found_motifs[idx]
    # df_found_motifs = df_found_motifs.loc[:,['motif_id', 'chr', 'start', 'end', 'cell_type', 'motif_id_without_orientation']]

    # motifs_list = list(to_reverse.keys())
    # motifs_list = [motif_id.split('_')[0] for motif_id in motifs_list]
    # motifs_list = list(set(motifs_list))

    # def find_pair_orientation(motif_id: str, row, df_found_motifs: pd.DataFrame):

    #     chromosome, x_start, x_end, y_start, y_end = row
    #     x_motif = df_found_motifs.loc[(df_found_motifs.loc[:,'chr']==chromosome)&(df_found_motifs.loc[:,'start']==x_start)&(df_found_motifs.loc[:,'end']==x_end)&(df_found_motifs.loc[:,'motif_id_without_orientation']==motif_id), 'motif_id']
    #     y_motif = df_found_motifs.loc[(df_found_motifs.loc[:,'chr']==chromosome)&(df_found_motifs.loc[:,'start']==y_start)&(df_found_motifs.loc[:,'end']==y_end)&(df_found_motifs.loc[:,'motif_id_without_orientation']==motif_id), 'motif_id']
    #     if len(x_motif) == 0 and len(y_motif) == 0:
    #         return 0
    #     elif len(x_motif) == 0 or len(y_motif) == 0:
    #         return 1
    #     else:
    #         x_motif_orientation = x_motif.to_string().split('_')[-1]
    #         y_motif_orientation = y_motif.to_string().split('_')[-1]
    #         if x_motif_orientation == y_motif_orientation:
    #             return 4
    #         elif x_motif_orientation == 'f' and y_motif_orientation == 'r':
    #             return 2
    #         elif x_motif_orientation == 'r' and y_motif_orientation == 'f':
    #             return 3

    # main_df_new = pl.from_pandas(main_df_new)
    # for motif_id in motifs_list:
    #     main_df_new = main_df_new.with_columns(main_df_new.select(['chr', 'x_start', 'x_end', 'y_start', 'y_end']).apply(
    #         lambda row: find_pair_orientation(motif_id, row, df_found_motifs),
    #         ).rename({'apply': f'{motif_id}_orientation'}))
    # main_df_new = main_df_new.to_pandas()

    # # Sum and join columns in main_df_new for the same motif on the same anchor (x or y)
    # for motif_id in motifs_list:
    #     to_drop = [f'x_{motif_id}_f', f'x_{motif_id}_r', f'y_{motif_id}_f', f'y_{motif_id}_r']
    #     motif_count = []
    #     for motif_anchor_orientation in to_drop:
    #         if motif_anchor_orientation in main_df_new.columns:
    #             motif_count.append(main_df_new[motif_anchor_orientation])
    #         else:
    #             motif_count.append(pd.Series([0]*len(main_df_new), dtype='int16'))

    #     main_df_new[f'x_{motif_id}_count'] = motif_count[0] + motif_count[1]
    #     main_df_new[f'y_{motif_id}_count'] = motif_count[2] + motif_count[3]
    #     main_df_new.drop(to_drop, axis=1, inplace=True)

    return main_df_new


def count_motifs(
    main_dfs_dict: dict, motifs_df: pd.DataFrame, organism: str = "human"
) -> Dict[str, pd.DataFrame]:
    """
    Counts motif occurances in both regions of each chromatin loop,
    for each dataframe from the main_dfs_dict dictionary.
    Args:
        main_dfs_dict: dictionary with cell types as keys and load functions of pandas DataFrames with chromatin loops as values.
        motifs_df: dictionary with cell types as keys and pandas DataFrame with motif occurrences found as values.
        organism: name of the organism.
    Returns:
        dictionary:
            keys: cell types
            values: pandas DataFrame with added columns of each motif counts in both regions of each chromatin loop
    """
    print("Adding motifs counts...")
    motifs_df[["chr", "cell_type"]] = motifs_df[["chr", "cell_type"]].astype("string")
    # if main_dfs_dict values are not DataFrames, load them
    if not isinstance(list(main_dfs_dict.values())[0], pd.DataFrame):
        main_dfs_dict = _dict_partitions(main_dfs_dict)

    motifs_groups = motifs_df.groupby("cell_type")

    for cell_type, motifs_df_group in motifs_groups:
        print(f"...for {cell_type} cell...")
        for main_name, main_df in main_dfs_dict.items():
            if main_name == cell_type:
                motifs_df_group = motifs_df_group.loc[
                    :, motifs_df_group.columns != "cell_type"
                ]
                main_df = _count_motifs_single_df(
                    main_df, motifs_df_group, cell_type, organism=organism
                )
                main_dfs_dict[main_name] = main_df
    print("Done!")

    return main_dfs_dict


def _remove_overlapping_single_df(df: pd.DataFrame) -> pd.DataFrame:
    """
    Remove negative examples that overlap with positive examples.
    Args:
        df: pandas DataFrame with positive and negative examples.
    Returns:
        pandas DataFrame with removed negative examples that overlap with positive examples.
    """
    df["id"] = [i for i in range(len(df))]
    df_1x = df.loc[df.loc[:, "label"] == 1, ["chr", "x_start", "x_end", "id"]]
    df_0x = df.loc[df.loc[:, "label"] == 0, ["chr", "x_start", "x_end", "id"]]
    df_1y = df.loc[df.loc[:, "label"] == 1, ["chr", "y_start", "y_end", "id"]]
    df_0y = df.loc[df.loc[:, "label"] == 0, ["chr", "y_start", "y_end", "id"]]
    # Find overlapping regions for x and y anchors separately
    x_overlap = _get_overlapping_regions(
        df_1x,
        df_0x,
        names=["chr1", "satrt1", "end1", "id_pos", "chr2", "start2", "end2", "id_neg"],
        wa=True,
        wb=True,
    )
    y_overlap = _get_overlapping_regions(
        df_1y,
        df_0y,
        names=["chr1", "satrt1", "end1", "id_pos", "chr2", "start2", "end2", "id_neg"],
        wa=True,
        wb=True,
    )
    # Find overlapping regions for x and y anchors together (rows with the same id_pos and id_neg)
    x = set(x_overlap.loc[:, ["id_pos", "id_neg"]].apply(tuple, axis=1))
    y = set(y_overlap.loc[:, ["id_pos", "id_neg"]].apply(tuple, axis=1))
    x_and_y = set([i[1] for i in x & y])
    # Remove overlapping negative examples
    df = df.loc[~df["id"].isin(x_and_y), :]
    df = df.loc[:, df.columns != "id"]

    assert len(df_1x) == len(
        df.loc[df.loc[:, "label"] == 1, :]
    ), "Something went wrong with removing overlapping negative examples"
    assert len(df_1y) == len(
        df.loc[df.loc[:, "label"] == 1, :]
    ), "Something went wrong with removing overlapping negative examples"

    return df


def _check_and_drop_nan_rows(df: pd.DataFrame) -> pd.DataFrame:
    """
    Checks if there are any rows with NaN values in the dataframe.
    Args:
        df: pandas DataFrame.
    Returns:
        pandas DataFrame without rows with NaN values.
    """
    print("\nChecking for NaN values in the dataframe...")
    if df.isnull().values.any():
        before = len(df)
        columns_with_na = df.columns[df.isna().any()].tolist()
        print(
            f"There are NaN values in the dataframe in columns: {columns_with_na}. Dropping rows with NaN values..."
        )
        df = df.dropna()
        after = len(df)
        print(f"Done! {before - after} rows were dropped.")
    else:
        print("There are no NaN values in the dataframe.")

    return df


def remove_overlapping(main_dfs_dict: dict):
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
    print("Removing overlapping negative examples...")
    # if main_dfs_dict values are not DataFrames, load them
    if not isinstance(list(main_dfs_dict.values())[0], pd.DataFrame):
        main_dfs_dict = _dict_partitions(main_dfs_dict)

    for main_name, main_df in main_dfs_dict.items():
        print(f"...for {main_name} cell...")
        len_before = len(main_df)
        main_df = _remove_overlapping_single_df(main_df)
        print(f"...{len_before - len(main_df)} examples were removed.")
        main_df = _check_and_drop_nan_rows(main_df)
        main_dfs_dict[main_name] = main_df
    print("Done!")

    return main_dfs_dict


def _set_cell_type(df: pd.DataFrame, cell_type: str) -> pd.DataFrame:
    """
    Set cell type in the dataframe.
    Args:
        df: pandas DataFrame.
        cell_type: name of the cell type.
    Returns:
        pandas DataFrame with set cell type.
    """
    df["cell_type"] = cell_type

    return df


def _drop_negs_with_same_anchors(df: pd.DataFrame) -> pd.DataFrame:
    df["chr;x"] = [df.iloc[i, 0] + ";" + str(df.iloc[i, 1]) for i in range(len(df))]
    df["chr;y"] = [df.iloc[i, 0] + ";" + str(df.iloc[i, 4]) for i in range(len(df))]
    df.reset_index(inplace=True)
    positive_x = set(df[df["label"] == 1]["chr;x"])
    negative_x = set(df[df["label"] == 0]["chr;x"])
    common_x = list(negative_x & positive_x)
    print("Common x number: " + str(len(common_x)))
    indexes_to_drop = df[(df["label"] == 0) & (df["chr;x"].isin(common_x))].index
    print("Shape before x drop: " + str(df.shape))
    df.drop(indexes_to_drop, axis=0, inplace=True)
    print("Shape after x drop: " + str(df.shape))

    positive_y = set(df[df["label"] == 1]["chr;y"])
    negative_y = set(df[df["label"] == 0]["chr;y"])
    common_y = list(negative_y & positive_y)
    print("Common y number: " + str(len(common_y)))
    indexes_to_drop = df[(df["label"] == 0) & (df["chr;y"].isin(common_y))].index
    print("Shape before y drop: " + str(df.shape))
    df.drop(indexes_to_drop, axis=0, inplace=True)
    print("Shape after y drop: " + str(df.shape))
    df.drop(["chr;x", "chr;y", "index"], axis=1, inplace=True)

    return df


def _random_negatives(df: pd.DataFrame):
    """
    Randomly select negative examples from other cell types.
    Args:
        df: pandas DataFrame.
        cell_type: name of the cell type.
    Returns:
        pandas DataFrame with randomly selected negative examples.
    """
    pos_num = len(df[df["label"] == 1])
    neg_num = len(df[df["label"] == 0])
    pos_indexes = df[df["label"] == 1].index
    neg_indexes = (
        df[df["label"] == 0].sample(frac=pos_num / neg_num, replace=False).index
    )

    df = df.loc[pos_indexes.union(neg_indexes), :]

    return df


def _random_pairing_negatives(
    df: pd.DataFrame, df_neg: pd.DataFrame, neg_pos_ratio: int, random_state: int = 42
):
    """
    Random pairing of negatives lying on the same chromosome from the df_neg dataframe.
    Filling with features from df dataframe.
    Args:
        df: pandas DataFrame with features.
        df_neg: pandas DataFrame with negative anchors.
        neg_pos_ratio: ratio of positive to negative examples.
        random_state: random state.
    """
    # random pairing of negatives lying on the same chromosome
    df_neg_pairs = df_neg[["chr", "center"]]
    df_neg_pairs = pl.from_pandas(df_neg_pairs)
    df_neg_pairs = df_neg_pairs.join(
        df_neg_pairs, how="outer", on="chr", suffix="_y"
    ).filter(pl.col("center") < pl.col("center_y"))
    df_neg_pairs = df_neg_pairs.rename({"center": "x", "center_y": "y"})
    df_neg_pairs = df_neg_pairs.to_pandas()
    df_neg_pairs = df_neg_pairs[df_neg_pairs["x"] < df_neg_pairs["y"]]
    # randomy sample negative examples WITHOUT REPEATS to get the desired neg_pos_retio
    positives_num = len(df[df["label"] == 1])
    if neg_pos_ratio:
        df_neg_pairs = df_neg_pairs.sample(
            n=round(int(neg_pos_ratio * positives_num)),
            random_state=random_state,
            replace=False,
        )
    # add features from df to df_neg_pairs
    df_neg = df[df["label"] == 0]
    df_neg = df_neg.drop_duplicates()
    df = df[df["label"] == 1]
    columns_to_add = list(df.columns)[9:]
    df_neg_pairs

    return df_neg_pairs


def concat_dfs_from_dict(
    main_dfs_dict: dict,
    cells_to_use: list = None,
    negtype=None,
    organism: str = "human",
) -> pd.DataFrame:
    """
    Concatenates dataframes from dictionary.
    Args:
        main_dfs_dict: dictionary with cell types as keys and load functions of pandas DataFrames with chromatin loops as values.
        cells_to_use: list of cell types to be used. If empty, all cell types from main_dfs_dict will be used.
        organism: name of the organism.
        negtype: type of negative examples to be used.
    Returns:
        pandas DataFrame with concatenated dataframes from dictionary.
    """
    cells_to_use = cells_to_use or []
    print("Concatenating dataframes from dictionary...")
    # check if all keys_to_concat are in main_dfs_dict
    assert set(cells_to_use) == set(
        main_dfs_dict.keys()
    ), "Something went wrong when filtering out the cell types to be used. Check data_preprocessing.yml file."

    for i, cell in enumerate(cells_to_use):
        print(cell)
        if i == 0:
            main_df = _set_cell_type(main_dfs_dict[cell], cell)
            if "distance" not in main_df.columns:
                main_df["distance"] = main_df["y"] - main_df["x"]
        else:
            df_to_concat = _set_cell_type(main_dfs_dict[cell], cell)
            if "distance" not in df_to_concat.columns:
                df_to_concat["distance"] = df_to_concat["y"] - df_to_concat["x"]
            main_df = pd.concat([main_df, df_to_concat], ignore_index=True)
        main_dfs_dict.pop(cell)

    main_df = _sort_df(main_df, "x_start", organism=organism)

    # print size of the dataframe
    print(f"Done! Shape of the dataframe: {main_df.shape}")

    return main_df
