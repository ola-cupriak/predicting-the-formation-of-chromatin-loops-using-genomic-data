import copy
import numpy as np
import math
import pyBigWig
from predicting_the_formation_of_chromatin_loops_using_genomic_data.utils import _dict_partitions
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



def _sort_df(df: pd.DataFrame, region_name: str) -> pd.DataFrame:
    """
    Sort the dataframe by chromosome and position.
    Args:
        df: pandas DataFrame.
        region_name: name of the column with regions.
    Returns:
        sorted pandas DataFrame.
    """
    df.loc[df['chr'] == 'X', 'chr'] = '100'
    df.loc[df['chr'] == 'Y', 'chr'] = '200'
    df['chr'] = df['chr'].astype(int)
    df = df.sort_values(by=['chr', region_name])
    df.loc[df['chr'] == 100, 'chr'] = 'X'
    df.loc[df['chr'] == 200, 'chr'] = 'Y'
    df['chr'] = df['chr'].astype('string')

    return df

    
def read_hic(partitioned_input: Dict[str, Callable[[], Any]], 
                cells2names: Dict[str, dict],
                dataset_name: str, r: int,
                cells_to_use: list=[]) -> Dict[str, pd.DataFrame]:
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
    dfs_dict = _dict_partitions(partitioned_input)
    cells2names_dataset_dict = cells2names[dataset_name]
    keys_dict = {".".join(key.split(".")[:-1]): key for key in cells2names_dataset_dict}
    assert set(cells_to_use).issubset(set(cells2names_dataset_dict.values())), f"Cell types: {set(cells_to_use)-set(cells2names_dataset_dict.values())} are not in the dataset. Please check data_preprocessing.yml file."
    new_dfs_dict = dict()
    for name, df in dfs_dict.items():
        cell_type = cells2names_dataset_dict[keys_dict[name]]
        if cells_to_use and cell_type not in cells_to_use:
            continue
        f1 = lambda x: x['chr1'].split("chr")[-1]
        f2 = lambda x: x['chr2'].split("chr")[-1]
        df["chr1"] = df.apply(f1, axis=1)
        df["chr2"] = df.apply(f2, axis=1)
        assert len(df[df['chr1']!=df['chr2']]) == 0
        df.rename(columns={'chr1': 'chr'}, inplace=True)
        for region in ['x', 'y']:
            df[f'{region}'] = (round((df[f'{region}1'] + df[f'{region}2'])/2)).astype(int)
            df[f'{region}_start'] = (df[f'{region}'] - r).astype(int)
            df[f'{region}_end'] = (df[f'{region}'] + r).astype(int)
        df['len_anchors'] = (df['x2']-df['x1'])+(df['y2']-df['y1']) # add length of both anchors
        df['len_loop'] = (df['y2']-df['x1']) # add length of loop
        df['cell_type'] = cell_type
        df = df[df.columns.intersection(['x', 'y', 'x_start', 'x_end', 'y_start', 'y_end', 'chr', 'cell_type'])]
        # sort by chr and x
        df = _sort_df(df, 'x')
        # set dtypes
        df = df.astype({'x': 'int32', 'y': 'int32', 'x_start': 'int32', 'x_end': 'int32', 
                        'y_start': 'int32', 'y_end': 'int32', 'chr': 'string', 'cell_type': 'string'})

        new_dfs_dict[cell_type] = df

    return new_dfs_dict


def read_peaks(partitioned_input: Dict[str, Callable[[], Any]],
                cells2names: Dict[str, dict],
                dataset_name: str,
                cells_to_use: list) -> Dict[str, pd.DataFrame]:
    """
    Load dataframes with DNase-seq/ChIP-seq peaks and modify the chromosome columns.
    Args:
        partitioned_input: dictionary with partition ids as keys and load functions as values.
        cells2names: dictionary, template: {'dataset_name': {'file_name': 'cell_type'}}
        dataset_name: name of the dataset to select from cells2names.
        cells_to_use: list of cell types to use.
    Returns:
        dictionary:
            keys: cell types 
            values: pandas DataFrames with DNase-seq/ChIP-seq peaks.
    """
    dfs_dict = _dict_partitions(partitioned_input)
    cells2names_dataset_dict = cells2names[dataset_name]
    assert set(cells_to_use).issubset(set(cells2names_dataset_dict.values())), f"Cell types: {set(cells_to_use)-set(cells2names_dataset_dict.values())} are not in the dataset. Please check data_preprocessing.yml file."
    keys_dict = {".".join(key.split(".")[:-1]): key for key in cells2names_dataset_dict}
    new_dfs_dict = dict()
    for name, df in dfs_dict.items():
        df = pl.from_pandas(df)
        cell_type = cells2names_dataset_dict[keys_dict[name]]
        if cells_to_use and cell_type not in cells_to_use:
            continue
        df = df.with_columns(df.select(['chr']).apply(
            lambda x: x[0].split('chr')[-1], return_dtype=pl.Utf8
            ).rename({'apply': 'chr'}))

        df = df.with_columns(pl.lit(cell_type).cast(pl.Utf8).alias('cell_type'))
        df = df.to_pandas()
        df = _sort_df(df, 'start')

        new_dfs_dict[cell_type] = df

    return new_dfs_dict


def read_bigWig(partitioned_input: Dict[str, Callable[[], Any]],
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
    input_dict = _dict_partitions(partitioned_input)
    cells2names_dataset_dict = cells2names[dataset_name]
    assert set(cells_to_use).issubset(set(cells2names_dataset_dict.values())), f"Cell types: {set(cells_to_use)-set(cells2names_dataset_dict.values())} are not in the dataset. Please check data_preprocessing.yml file."
    keys_dict = {".".join(key.split(".")[:-1]): key for key in cells2names_dataset_dict}
    bigWig_data_dict = dict()
    for name, bigWig_path in input_dict.items():
        cell_type = cells2names_dataset_dict[keys_dict[name]]
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


def _get_overlapping_regions(df1: pd.DataFrame, df2: pd.DataFrame, names: list, count: bool = False, 
                            wa: bool = False, wb: bool = False) -> pd.DataFrame:
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


# def _remove_neg_overlapping_pos(pos_df, neg_df):
#     pos_df['id'] = [i for i in range(len(pos_df))]
#     neg_df['id'] = [i for i in range(len(neg_df))]
#     # Find overlapping regions for x and y anchors separately
#     x_overlap = _get_overlapping_regions(pos_df[['chr', 'x_start', 'x_end', 'id']], neg_df[['chr', 'x_start', 'x_end', 'id']], names=['chr1', 'satrt1', 'end1', 'id_pos', 'chr2', 'start2', 'end2', 'id_neg'], wa=True, wb=True)
#     y_overlap = _get_overlapping_regions(pos_df[['chr', 'y_start', 'y_end', 'id']], neg_df[['chr', 'y_start', 'y_end', 'id']], names=['chr1', 'satrt1', 'end1', 'id_pos', 'chr2', 'start2', 'end2', 'id_neg'], wa=True, wb=True)
#     # Find overlapping regions for x and y anchors together (rows with the same id_pos and id_neg)
#     x = set(x_overlap.loc[:,['id_pos', 'id_neg']].apply(tuple, axis=1))
#     y = set(y_overlap.loc[:,['id_pos', 'id_neg']].apply(tuple, axis=1))
#     x_and_y = set([i[1] for i in x & y])
#     # Remove overlapping negative examples
#     neg_df = neg_df.loc[~neg_df['id'].isin(x_and_y),:]
#     neg_df = neg_df.loc[:, neg_df.columns != 'id']

#     return neg_df


def _create_pairs_each_with_each_single_df(df_pos, df_open_chromatin, cell_type, r, neg_pos_ratio, random_state) -> pd.DataFrame:
    """
    """
    df_open_chromatin['center'] = (df_open_chromatin['start'] + df_open_chromatin['end']) // 2
    df_open_chromatin = df_open_chromatin[['chr', 'center']] 
    chromosomes = pd.unique(df_pos['chr'])
    chromosomes = {k: round(len(df_pos[df_pos['chr'] == k])*neg_pos_ratio) for k in chromosomes}
    df_neg = pd.DataFrame()

    for ch in chromosomes.keys():
        df_chr = df_open_chromatin[df_open_chromatin['chr'] == ch]
        df_chr = pl.from_pandas(df_chr)
        if len(df_chr) > 1000:
            df_chr = df_chr.sample(n=1000, with_replacement=False, seed=random_state, shuffle=True)
        df_chr = df_chr.join(df_chr, how='outer', on='chr', suffix='_y')
        # rename columns in polars
        df_chr = df_chr.rename({'center': 'x', 'center_y': 'y'})
        # filter columns in polars
        df_chr = df_chr.filter(pl.col('x') < pl.col('y'))
        df_chr = df_chr.filter(pl.col('x') + r < pl.col('y') - r)
        # sample with polars
        df_chr = df_chr.sample(n=chromosomes[ch], with_replacement=False, seed=random_state, shuffle=True)
        # convert back to pandas
        df_chr = df_chr.to_pandas()
        df_neg = pd.concat([df_neg, df_chr])
    
    df_neg['cell_type'] = cell_type
    df_neg['label'] = 0
    df_neg['x_start'] = df_neg['x'] - r
    df_neg['x_end'] = df_neg['x'] + r
    df_neg['y_start'] = df_neg['y'] - r
    df_neg['y_end'] = df_neg['y'] + r

    df_pairs = pd.concat([df_pos, df_neg])
    df_pairs.reset_index(drop=True, inplace=True)

    df_pairs = df_pairs.astype({'x': 'int32', 'y': 'int32', 'x_start': 'int32', 'x_end': 'int32', 
                                'y_start': 'int32', 'y_end': 'int32', 'label': 'int16'})
    
    return df_pairs


def _x_and_y_anchors2one_col(df: pd.DataFrame, 
                             x_cols2use: list=['chr', 'x_start', 'x_end'], 
                             y_cols2use: list=['chr', 'y_start', 'y_end'], 
                             new_cols: list=['chr', 'start', 'end'], 
                             sortby:str='start') -> pd.DataFrame:
    """
    Combines the columns describing x regions and the columns describing y regions 
    into one set of columns describing all regions.
    (columns: x_chr, x_start, x_end, y_chr, y_start, y_end -> columns: chr, start, end)
    Args:
        df: DataFrame with columns describing x regions and columns describing y regions.
    Returns:
        pandas DataFrame with one set of columns describing all regions.
    """
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
    anchors_df = _sort_df(anchors_df, sortby)
    # reset index
    anchors_df = anchors_df.reset_index(drop=True)
    anchors_df = anchors_df[new_cols]

    return anchors_df


def _create_new_anchors_pairs(df: pd.DataFrame):
    """
    Merge all values in one column (x) with all values in second column (y)
    """
    df_new_pairs = _x_and_y_anchors2one_col(df, ['chr', 'x'], ['chr', 'y'], ['chr', 'anchor'], 'anchor')
    df_new_pairs = pl.from_pandas(df_new_pairs)
    df_new_pairs = df_new_pairs.join(df_new_pairs, how='outer', on='chr', suffix='_y').filter(pl.col('anchor') < pl.col('anchor_y'))
    df_new_pairs = df_new_pairs.rename({'anchor': 'x', 'anchor_y': 'y'})
    df_new_pairs = df_new_pairs.to_pandas()
    
    return df_new_pairs
    

def _remove_duplicated_pairs(df_pos: pd.DataFrame, df_neg: pd.DataFrame) -> pd.DataFrame:
    """
    Remove pairs from neg_df that are present in pos_df
    Args:
        df_pos: pandas DataFrame with positive examples.
        df_neg: pandas DataFrame with negative examples.
    Returns:
        pandas DataFrame with negative examples that are not present in positive examples.
    """
    df_pos = df_pos[['chr', 'x', 'y']]
    df_neg = df_neg[['chr', 'x', 'y']]
    df_neg = df_neg.merge(df_pos, how='left', on=['chr', 'x', 'y'], indicator=True)
    df_neg = df_neg[df_neg['_merge'] == 'left_only']
    df_neg.drop(columns='_merge', inplace=True)

    return df_neg


def _get_negatives_by_new_anchors_pairing(df: pd.DataFrame, cell_type: str, r: int, neg_pos_ratio: float, random_state: int, df_len: int):
    df_neg = _create_new_anchors_pairs(df)
    df_neg = _remove_duplicated_pairs(df, df_neg)
    df_neg['x_start'] = df_neg['x'] - r
    df_neg['x_end'] = df_neg['x'] + r
    df_neg['y_start'] = df_neg['y'] - r
    df_neg['y_end'] = df_neg['y'] + r
    df_neg = df_neg[df_neg['x_end'] < df_neg['y_start']]
    df_neg['cell_type'] = cell_type
    df_neg['label'] = 0
    df_neg = df_neg[['chr', 'x', 'x_start', 'x_end', 'y', 'y_start', 'y_end', 'cell_type', 'label']]
    # randomy sample negative examples WITHOUT REPEATS to get the desired neg_pos_retio
    df_neg = df_neg.sample(n=round(int(neg_pos_ratio * df_len)), random_state=random_state, replace=False)
    # merge positive and negative examples
    df = pd.concat([df, df_neg], axis=0)
    # sort by chr and region
    df = _sort_df(df, 'x_start')
    df = df.reset_index(drop=True)
    df = df.astype({'x': 'int32', 'y': 'int32', 'x_start': 'int32', 'x_end': 'int32', 
                    'y_start': 'int32', 'y_end': 'int32', 'label': 'int16'})

    return df


def add_labels(dfs_dict: Dict[str, pd.DataFrame], mtype: str, peaks_dict: Dict[str, pd.DataFrame], r, neg_pos_ratio: float, random_state: int) -> None:
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
    assert mtype in ['anchors_from_other_cell_types', 
                    'new_pairs_of_anchors', 
                    'open_chromatin_regions',
                    'new_pairs_and_open_chromatin'], f"Negative sampling type should be either 'anchors_from_other_cell_types', 'new_pairs_of_anchors' or 'open_chromatin_regions, but got {type}."
    if mtype == 'anchors_from_other_cell_types':
        df = _concat_dfs(dfs_dict)
        df = _sort_df(df, 'x')
    
    for name, cell_df in dfs_dict.items():
        print('Creating negative examples for cell type', name, '...')

        if mtype == 'anchors_from_other_cell_types':
            f = lambda x: 1 if x['cell_type'] == name else 0
            df['label'] = df.apply(f, axis=1)
            df['label'] = df['label'].astype('int16')
            dfs_dict[name] = copy.deepcopy(df)
        elif mtype == 'new_pairs_of_anchors':
            cell_df['label'] = 1
            df_len = len(cell_df)
            dfs_dict[name] = _get_negatives_by_new_anchors_pairing(cell_df, name, r, neg_pos_ratio, random_state, df_len)
        elif mtype == 'open_chromatin_regions':
            cell_df['label'] = 1
            dfs_dict[name] = _create_pairs_each_with_each_single_df(cell_df, peaks_dict[name], name, r, neg_pos_ratio, random_state)
        else:
            cell_df['label'] = 1
            df_len = len(cell_df)
            with_neg = _create_pairs_each_with_each_single_df(cell_df, peaks_dict[name], name, r, neg_pos_ratio/2, random_state)
            with_neg = _get_negatives_by_new_anchors_pairing(with_neg, name, r, neg_pos_ratio/2, random_state, df_len)
            dfs_dict[name] = with_neg
    
    return dfs_dict


def _find_the_closest_peak(main_df: pd.DataFrame, peak_df: pd.DataFrame) -> pd.DataFrame:
    """
    Find the distance of the closest peak from the centre of each region.
    Args:
        main_df: pandas DataFrame with regions.
        peak_df: pandas DataFrame with peaks.
    Returns:
        pandas DataFrame with the distances of the closest peaks to the centre of each region.
    """ 
    # Remove peaks on chromosomes that are not in the main_df
    chromosomes = main_df['chr'].unique()
    peak_df = peak_df.loc[peak_df['chr'].isin(chromosomes),:]

    assert len(main_df.loc[main_df.iloc[:,1] != main_df.iloc[:,2],:]) == 0, 'The region columns are not equal'
    bed1 = pybedtools.BedTool.from_dataframe(main_df)
    bed2 = pybedtools.BedTool.from_dataframe(peak_df)

    closest = bed1.closest(bed2, d=True, t='first')
    closest = pd.read_table(closest.fn, names=['chr', 'centre', 'centre_dup', 'second_reg', 'cell_type',
                                            'chr_peak', 'start_peak', 'end_peak', 'cell_type_peak', 'distance'])
    closest = closest.loc[:,['chr', 'centre', 'second_reg', 'distance', 'cell_type']]

    return closest


def _count_peaks_single_df(main_df: pd.DataFrame, peaks_df: pd.DataFrame, experiment: str) -> pd.DataFrame:
    """
    Count the number of experiment peaks in both regions of each chromatin loop.
    Args:
        main_df: pandas DataFrame with chromatin loops.
        peaks_df: pandas DataFrame with experiment peaks.
        experiment: name of the experiment.
    Returns:
        pandas DataFrame with added columns of numbers of experiment peaks in both regions of each loop
    """
    for region in ['x', 'y']:
        second_reg = 'y' if region == 'x' else 'x'
        peak_counts = main_df.loc[:,['chr', f'{region}_start', f'{region}_end', f'{second_reg}', 'cell_type']]
        peak_counts = _get_overlapping_regions(peak_counts, peaks_df, names=['chr', 'start', 'end', 'second_reg', 'cell_type', 'count'], count=True)
        peak_counts.rename(columns={'count': f'{region}_{experiment}_counts',
                                    'start': f'{region}_start',
                                    'end': f'{region}_end',
                                    'second_reg': f'{second_reg}'}, inplace=True)

        assert len(peak_counts) == len(main_df), 'Length of the main_df and peak_counts are not equal'

        main_len_before = len(main_df)
        main_df = main_df.merge(peak_counts)

        assert len(main_df) == main_len_before, 'Length of the main_df changed after merging'

    to_change_dtype = [f'x_{experiment}_counts', f'y_{experiment}_counts']
    main_df.loc[:,to_change_dtype] = main_df.loc[:,to_change_dtype].astype('int16')

    return main_df


def _find_the_closest_peaks_single_df(main_df: pd.DataFrame, peaks_df: pd.DataFrame, experiment: str) -> pd.DataFrame:
    """
    Find the distance of the closest peak of the experiment from the center of each anchor from the chromatin loop.
    Args:
        main_df: pandas DataFrame with chromatin loops.
        peaks_df: pandas DataFrame with experiment peaks.
        experiment: name of the experiment.
    Returns:
        pandas DataFrame with added columns of distances of the closest peak of the experiment 
        from the center of each anchor from the chromatin loop.
    """
    for region in ['y', 'x']:
        second_reg = 'y' if region == 'x' else 'x'
        distances = main_df.loc[:,['chr', f'{region}', f'{second_reg}', 'cell_type']]
        distances = _sort_df(distances, region)
        distances = distances.loc[:,['chr', f'{region}', f'{region}', f'{second_reg}', 'cell_type']]
        distances = _find_the_closest_peak(distances, peaks_df)
        distances.rename(columns={'distance': f'{region}_{experiment}_distance',
                                    'centre': f'{region}',
                                    'second_reg': f'{second_reg}'
                                    }, inplace=True)

        assert len(distances) == len(main_df), 'Length of the main_df and distances are not equal'

        main_len_before = len(main_df)
        main_df = main_df.merge(distances)

        assert len(main_df) == main_len_before, 'Length of the main_df changed after merging'

    to_change_dtype = [f'x_{experiment}_distance', f'y_{experiment}_distance']
    main_df.loc[:,to_change_dtype] = main_df.loc[:,to_change_dtype].astype('int16')

    return main_df


def count_peaks_and_distances(main_dfs_dict: Dict[str, Callable[[], Any]], peaks_dfs_dict: Dict[str, pd.DataFrame], 
                experiment: str) -> pd.DataFrame:
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
    print(f'Adding peaks counts and distances from anchor centers for {experiment}...')
    # if main_dfs_dict values are not DataFrames, load them
    if not isinstance(list(main_dfs_dict.values())[0], pd.DataFrame):
        main_dfs_dict = _dict_partitions(main_dfs_dict)
    for peaks_name, peaks_df in peaks_dfs_dict.items():
        print(f'...for {peaks_name} cell...')
        for main_name, main_df in main_dfs_dict.items():
            if main_name == peaks_name:
                main_df = _count_peaks_single_df(main_df, peaks_df, experiment)
                main_df = _find_the_closest_peaks_single_df(main_df, peaks_df, experiment)
                main_dfs_dict[main_name] = main_df
    print('Done!')

    return main_dfs_dict


def get_triangle_kernel(kerlen: int) -> list:
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


def calculate_weighted_mean(distribution: list):
    """
    Calculate weighted mean of the distribution.
    Args:
        distribution: list with values of the distribution.
    Returns:
        weighted mean of the distribution.
    """
    weights = get_triangle_kernel(len(distribution))
    return sum([distribution[i]*weights[i] for i in range(len(distribution))]) / sum(weights)


def _add_bigWig_data_single_df(main_df: pd.DataFrame, bigWig_path, experiment: str) -> pd.DataFrame:
    """
    POLARS
    Count statistics (weighted mean, arithmetic mean, minimum and maximum) 
    of the bigWig data in both regions of each chromatin loop.
    Args:
        main_df: pandas DataFrame with chromatin loops.
        bigWig_obj: pyBigWig object with experiment peaks
    Returns:
        pandas DataFrame with added columns of bigWig data statistics in both regions of each loop
    """
    bigWig_obj = pyBigWig.open(bigWig_path)
    regions = ['x', 'y']
    main_df = pl.from_pandas(main_df)
    for region in regions:
        main_df = main_df.with_columns(main_df.select(['chr', f'{region}_start', f'{region}_end']).apply(
            lambda x: calculate_weighted_mean(bigWig_obj.values('chr'+x[0], x[1], x[2])), return_dtype=pl.Float32
            ).rename({'apply':f'{region}_{experiment}_weighted_mean'}))
        
        main_df = main_df.with_columns(main_df.select(['chr', f'{region}_start', f'{region}_end']).apply(
            lambda x: bigWig_obj.stats('chr'+x[0], x[1], x[2], type='mean')[0], return_dtype=pl.Float32
            ).rename({'apply':f'{region}_{experiment}_mean'}))
    
        main_df = main_df.with_columns(main_df.select(['chr', f'{region}_start', f'{region}_end']).apply(
            lambda x: bigWig_obj.stats('chr'+x[0], x[1], x[2], type='max')[0], return_dtype=pl.Float32
            ).rename({'apply':f'{region}_{experiment}_max'}))
        
        main_df = main_df.with_columns(main_df.select(['chr', f'{region}_start', f'{region}_end']).apply(
            lambda x: bigWig_obj.stats('chr'+x[0], x[1], x[2], type='min')[0], return_dtype=pl.Float32
            ).rename({'apply':f'{region}_{experiment}_min'}))
        
    main_df = main_df.to_pandas()
    return main_df


def add_bigWig_data(main_dfs_dict: Dict[str, Callable[[], Any]],
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
    print(f'Adding {experiment} data...')
    # if main_dfs_dict values are not DataFrames, load them
    if not isinstance(list(main_dfs_dict.values())[0], pd.DataFrame):
        main_dfs_dict = _dict_partitions(main_dfs_dict)
    for bigWig_name, bigWig_path in bigWig_data_dict.items():
        print(f'...for {bigWig_name} cell...')
        for main_name, main_df in main_dfs_dict.items():
            if main_name == bigWig_name:
                main_df = _add_bigWig_data_single_df(main_df, bigWig_path, experiment)
                main_dfs_dict[main_name] = main_df
    print('Done!')

    return main_dfs_dict


def all_anchors2one_df(dfs_dict: Dict[str, pd.DataFrame]) -> pd.DataFrame:
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
    # if dfs_dict values are not DataFrames, load them
    if not isinstance(list(dfs_dict.values())[0], pd.DataFrame):
        dfs_dict = _dict_partitions(dfs_dict)

    df_with_2_regions = _concat_dfs(dfs_dict)
    anchors_df = _x_and_y_anchors2one_col(df_with_2_regions)
    
    return anchors_df


def all_peaks2one_df(peaks_dict: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    """
    Concatenates all peaks DataFrames into one DataFrame.
    Args:
        peaks_dict: dictionary with cell types as keys and pandas DataFrames as values.
    Returns:
        pandas DataFrame with all peaks.
    """
    df = _concat_dfs(peaks_dict)
    # sort by chr and region
    df = _sort_df(df, 'start')
    df = df.reset_index(drop=True)

    return df


def get_overlaps_with_names(anchors_df: pd.DataFrame, peaks_df: pd.DataFrame) -> pd.DataFrame:
    """
    Get overlapping regions between anchors and peaks.
    Args:
        anchors_df: pandas DataFrame with anchors coordinates.
        peaks_df: pandas DataFrame with peaks coordinates.
    Returns:
        pandas DataFrame with overlapping regions coordinates and names.
    """
    intersection_wa_wb = _get_overlapping_regions(anchors_df, peaks_df, names=['anchor_chr', 'anchor_start', 'anchor_end', 'peak_chr', 
                                                                              'peak_start', 'peak_end', 'cell_type'], wa=True, wb=True)
    intersection = _get_overlapping_regions(anchors_df, peaks_df, names=['chr', 'start', 'end'])

    joined_intersection = pd.merge(intersection_wa_wb, intersection, left_index=True, right_index=True)
    joined_intersection['name'] = joined_intersection.apply(lambda x: f"{x['chr']}:{x['anchor_start']}:{x['anchor_end']}:{x['cell_type']}", axis=1)

    return joined_intersection[['chr', 'start', 'end', 'name']]


def getfasta_bedfile(df: pd.DataFrame, path_simp_genome: str, path_to_save: str) -> str:
    """
    Cut sequences from chromosomes using BEDTools for coordinates from the pandas DataFrame.
    Args:
        df: pandas DataFrame with coordinates.
        path_simp_genome: path to fasta file with chromosomes.
        path_to_save: path to save fasta file.
    Returns:
        path_to_save: path to fasta file.
    """
    fasta = pybedtools.BedTool(path_simp_genome)
    df_to_search = df
    df_to_search['chr'] = 'chr' + df_to_search['chr'].astype(str)
    bed = pybedtools.BedTool.from_dataframe(df_to_search)
    fasta_bed = bed.sequence(fi=fasta, nameOnly=True)

    records = []
    for record in SeqIO.parse(fasta_bed.seqfn, 'fasta'):
        records.append(record)
    
    SeqIO.write(records, path_to_save, "fasta")

    return path_to_save


def _modify_output(line: str, i: int):
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
        line = line.replace('sequence name', 'chr\tanchor_start\tanchor_end\tcell_type')
        line = line.replace('#pattern name', 'motif_id')
    else:
        # Get the information you need 
        elements = line.split('\t')
        motif_id, name, _, _, strand, _, _, _, _ = elements
        # Save motif_id and motif_alt_id to the dictionary
        # Change name format to tsv format
        name = name.replace(':', '\t')
        # Add strand information to the motif_id
        if strand == '+':
            motif_id += '_f'
        else:
            motif_id += '_r'
        # Update the list of lines
        elements[0] = motif_id
        elements[1] = name
        line = '\t'.join(elements)
    
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
    path_for_meme = path_motifs.replace('.txt', '.meme')
    subprocess.run(f'jaspar2meme -bundle {path_motifs} > {path_for_meme}', shell=True)
    # Find motif occurances in fasta sequences
    start = time.time()
    print('Finding motifs...')
    subprocess.run(f'fimo --text {path_for_meme} {path_fasta} > data/temp/temp.csv', shell=True)
    print(f'Finding motifs took {time.time() - start} seconds')
    # Read output
    lines = []
    n_lines = 0
    with open('data/temp/temp.csv', 'r') as f:
        for i, line in tqdm(enumerate(f)):
            lines.append(_modify_output(line, i))
            n_lines = i
            if i % 1000000 == 0:
                print(i)
                # Save modified output to temporary file
                with open('data/temp/temp2.csv', 'a') as f2:
                    f2.writelines(lines)
                    del lines
                    lines = []
        with open('data/temp/temp2.csv', 'a') as f2:
            f2.writelines(lines)
            del lines
    # Read temporary file as pandas DataFrame
    dtypes = {'chr': "string", 'start': "int32", 'end': "int32", 'motif_id': "string", 'cell_type': "string"}
    df = pd.read_csv('data/temp/temp2.csv', sep='\t', dtype=dtypes, usecols=list(dtypes.keys()))
    subprocess.run('rm data/temp/temp.csv', shell=True)
    subprocess.run('rm data/temp/temp2.csv', shell=True)
    df = df[['chr', 'start', 'end', 'motif_id', 'cell_type']]
    len_before = len(df)
    # Count motif occurences
    df = pd.DataFrame(df.groupby(['chr', 'start', 'end', 'cell_type', 'motif_id'], observed=True).size(), columns=['count'])
    df['count'] = df['count'].astype('int16')
    df = df.unstack('motif_id', fill_value=0)
    df.columns = ['_'.join(x) for x in df.columns if x[0]=='count']
    df.reset_index(inplace=True)
    df.rename(columns={name: name.replace('count_', '') for name in df.columns}, inplace=True)
    
    assert sum(df.iloc[:, 4:].sum(axis=1)) == len_before, 'Something went wrong with counting motifs'
    
    return df
    
    # DIFFERENT APPROACH - NOT WORKING
    # dtypes = {'chr': "string", 'anchor_start': "int32", 'anchor_end': "int32", 'motif_id': "string", 'cell_type': "string"}
    # df = pd.read_csv('data/temp/temp2.csv', sep='\t', dtype=dtypes, usecols=list(dtypes.keys()))
    # #####subprocess.run('rm data/temp/temp.csv', shell=True)
    # df = df[['chr', 'anchor_start', 'anchor_end', 'motif_id', 'cell_type']]
    # len_before = len(df)
    # # Count motif occurences
    # df.rename(columns={'anchor_start': 'start', 'anchor_end': 'end'}, inplace=True)
    # motif_cols = list(pd.unique(df['motif_id']))
    # df = pl.from_pandas(df)
    # all_columns = ['chr', 'start', 'end', 'cell_type'] + motif_cols
    # # save to parquet
    # for motif_id in tqdm(motif_cols):
    #     df.with_columns([df['motif_id'] == motif_id]).write_parquet(f'data/temp/temp_{motif_id}.parquet')


    # df.to_parquet('data/temp/temp2.parquet', engine='pyarrow')
    # del df

    # df = pd.DataFrame(columns=all_columns)

    # parquet_file = pq.ParquetFile('data/temp/temp2.parquet')
    # for batch in parquet_file.iter_batches(batch_size=1000000):
    #     batch_df = batch.to_pandas()
    #     batch_df = pd.DataFrame(batch_df.groupby(['chr', 'start', 'end', 'cell_type', 'motif_id'], observed=True).size(), columns=['count'])
    #     batch_df['count'] = batch_df['count'].astype('int16')
    #     batch_df = batch_df.unstack('motif_id', fill_value=0)
    #     batch_df.columns = ['_'.join(x) for x in batch_df.columns if x[0]=='count']
    #     batch_df.reset_index(inplace=True)
    #     batch_df.rename(columns={name: name.replace('count_', '') for name in batch_df.columns}, inplace=True)
    #     # add non existing columns with zeros
    #     for col in all_columns:
    #         if col not in batch_df.columns:
    #             batch_df[col] = 0
    #     df = df.append(batch_df, ignore_index=True)
    
    # assert sum(df.iloc[:, 4:].sum(axis=1)) == len_before, 'Something went wrong with counting motifs'
    
    # return df


def _reverse_names(colnames: list) :
    """
    Create dictionary with swapped suffixes.
    Columns with suffix '_f' will be swapped to columns with suffix '_r' and vice versa.
    Args:
        colnames: list of column names.
    Returns:
        dictionary with swapped suffixes.
    """
    to_reverse_f = [name for name in colnames if '_f' in name]
    to_reverse_r = [name for name in colnames if '_r' in name]
    to_reverse_f = {name: name.replace('_f', '_r') for name in to_reverse_f}
    to_reverse_r = {name: name.replace('_r', '_f') for name in to_reverse_r}

    to_reverse = {**to_reverse_f, **to_reverse_r}

    return to_reverse


def _count_motifs_single_df(main_df: pd.DataFrame, motifs_df: pd.DataFrame) -> pd.DataFrame:
    """
    Counts motif occurances in both regions of each chromatin loop.
    Args:
        main_df: pandas DataFrame with chromatin loops.
        motifs_df: pandas DataFrame with motif occurrences found.
    Returns:
        pandas DataFrame with added columns of each motif counts in both regions of each chromatin loop.
    """
    # change dtypes in main_df
    main_df[['chr', 'cell_type']] = main_df[['chr', 'cell_type']].astype('string')

    # add id column
    main_df['id'] = [i for i in range(len(main_df))]
    main_df_new = copy.deepcopy(main_df)
    columns_no = len(main_df_new.columns)

    # reverse _f/_r suffix in names of columns with motifs
    motifs_colnames = list(motifs_df.columns)
    motifs_colnames.remove('chr')
    to_reverse = _reverse_names(motifs_colnames)

    for region in ['x', 'y']:
        if region == 'y':
            motifs_df = motifs_df.rename(columns=to_reverse)
        motifs_df_renamed = motifs_df.rename(columns={name: f'{region}_{name}' for name in motifs_colnames})
        main_df_new = main_df_new.merge(motifs_df_renamed, on=['chr', f'{region}_start', f'{region}_end'], how='inner')
    # add rows with missing ids 
    motifs_list = list(main_df_new.columns)[columns_no:]
    to_add = main_df[~main_df['id'].isin(main_df_new['id'])]
    for m in motifs_list:
        to_add.loc[:, m] = pd.Series([0]*len(main_df), dtype='int16')
    main_df_new = pd.concat([main_df_new, to_add], ignore_index=True)
    main_df_new = _sort_df(main_df_new, 'x_start')
    main_df_new = main_df_new.loc[:, main_df_new.columns != 'id']

    return main_df_new


def count_motifs(main_dfs_dict: dict, motifs_df: pd.DataFrame):
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
    print('Adding motifs counts...')
    motifs_df[['chr', 'cell_type']] = motifs_df[['chr', 'cell_type']].astype('string')
    # if main_dfs_dict values are not DataFrames, load them
    if not isinstance(list(main_dfs_dict.values())[0], pd.DataFrame):
        main_dfs_dict = _dict_partitions(main_dfs_dict)

    motifs_groups = motifs_df.groupby('cell_type')

    for cell_type, motifs_df_group in motifs_groups:
        print(f'...for {cell_type} cell...')
        for main_name, main_df in main_dfs_dict.items():
            if main_name == cell_type:
                motifs_df_group = motifs_df_group.loc[:, motifs_df_group.columns != 'cell_type']
                main_df = _count_motifs_single_df(main_df, motifs_df_group)
                main_dfs_dict[main_name] = main_df
    print('Done!')

    return main_dfs_dict


def _remove_overlapping_single_df(df: pd.DataFrame) -> pd.DataFrame:
    """
    Remove negative examples that overlap with positive examples.
    Args:
        df: pandas DataFrame with positive and negative examples.
    Returns:
        pandas DataFrame with removed negative examples that overlap with positive examples.
    """
    df['id'] = [i for i in range(len(df))]
    df_1x = df.loc[df.loc[:,'label']==1,['chr', 'x_start', 'x_end', 'id']]
    df_0x = df.loc[df.loc[:,'label']==0,['chr', 'x_start', 'x_end', 'id']]
    df_1y = df.loc[df.loc[:,'label']==1,['chr', 'y_start', 'y_end', 'id']]
    df_0y = df.loc[df.loc[:,'label']==0,['chr', 'y_start', 'y_end', 'id']]
    # Find overlapping regions for x and y anchors separately
    x_overlap = _get_overlapping_regions(df_1x, df_0x, names=['chr1', 'satrt1', 'end1', 'id_pos', 'chr2', 'start2', 'end2', 'id_neg'], wa=True, wb=True)
    y_overlap = _get_overlapping_regions(df_1y, df_0y, names=['chr1', 'satrt1', 'end1', 'id_pos', 'chr2', 'start2', 'end2', 'id_neg'], wa=True, wb=True)
    # Find overlapping regions for x and y anchors together (rows with the same id_pos and id_neg)
    x = set(x_overlap.loc[:,['id_pos', 'id_neg']].apply(tuple, axis=1))
    y = set(y_overlap.loc[:,['id_pos', 'id_neg']].apply(tuple, axis=1))
    x_and_y = set([i[1] for i in x & y])
    # Remove overlapping negative examples
    df = df.loc[~df['id'].isin(x_and_y),:]
    df = df.loc[:, df.columns != 'id']

    assert len(df_1x) == len(df.loc[df.loc[:,'label']==1,:]), 'Something went wrong with removing overlapping negative examples'
    assert len(df_1y) == len(df.loc[df.loc[:,'label']==1,:]), 'Something went wrong with removing overlapping negative examples'
    
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
    print('Removing overlapping negative examples...')
    # if main_dfs_dict values are not DataFrames, load them
    if not isinstance(list(main_dfs_dict.values())[0], pd.DataFrame):
        main_dfs_dict = _dict_partitions(main_dfs_dict)

    for main_name, main_df in main_dfs_dict.items():
        print(f'...for {main_name} cell...')
        len_before = len(main_df)
        main_df = _remove_overlapping_single_df(main_df)
        print(f'...{len_before - len(main_df)} examples were removed.')
        main_dfs_dict[main_name] = main_df
    print('Done!')

    return main_dfs_dict


def concat_dfs_from_dict(main_dfs_dict: dict, cells_to_use: list=[]) -> pd.DataFrame:
    '''
    Concatenates dataframes from dictionary.
    Args:   
        main_dfs_dict: dictionary with cell types as keys and load functions of pandas DataFrames with chromatin loops as values.
        cells_to_use: list of cell types to be used. If empty, all cell types from main_dfs_dict will be used.
    Returns:
        pandas DataFrame with concatenated dataframes from dictionary.
    '''
    print('Concatenating dataframes from dictionary...')
    # check if all keys_to_concat are in main_dfs_dict
    if cells_to_use:
        assert set(cells_to_use) == set(main_dfs_dict.keys()), 'Something went wrong when filtering out the cell types to be used. Check data_preprocessing.yml file.'
    else:
        cells_to_use = list(cells_to_use.keys())

    expected_len = sum([len(main_dfs_dict[key]) for key in cells_to_use])

    for i in range(len(cells_to_use)):
        if i == 0:
            main_df = main_dfs_dict[cells_to_use[i]]
        else:
            main_df = pd.concat([main_df, main_dfs_dict[cells_to_use[i]]], ignore_index=True)
        main_dfs_dict.pop(cells_to_use[i])
            
    assert len(main_df) == expected_len, 'Something went wrong with concatenating dataframes from dictionary - expected length is not equal to actual length.'

    main_df = _sort_df(main_df, 'x_start')

    print(f'Done!')

    return main_df