import pyBigWig
from predicting_the_formation_of_chromatin_loops_using_genomic_data.utils import _dict_partitions
from typing import Any, Callable, Dict
import pandas as pd
import copy
from tqdm import tqdm

def read_hic(partitioned_input: Dict[str, Callable[[], Any]], 
                cells2names: Dict[str, dict],
                dataset_name: str) -> Dict[str, pd.DataFrame]:
    """
    Load and modify the dataframes with chromatin loops anotations.
    Args:
        partitioned_input: dictionary with partition ids as keys and load functions as values.
        cells2names: dictionary, template: {'dataset_name': {'file_name': 'cell_type'}}
        dataset_name: name of the dataset to select from cells2names.
    Returns:
        dictionary with partition ids as keys and pandas DataFrames as values.
    """
    dfs_dict = _dict_partitions(partitioned_input)
    cells2names_dataset_dict = cells2names[dataset_name]
    keys_dict = {".".join(key.split(".")[:-1]): key for key in cells2names_dataset_dict}
    new_dfs_dict = dict()
    for name, df in dfs_dict.items():
        cell_type = cells2names_dataset_dict[keys_dict[name]]
        f1 = lambda x: x['chr1'].split("chr")[-1]
        f2 = lambda x: x['chr2'].split("chr")[-1]
        df["chr1"] = df.apply(f1, axis=1)
        df["chr2"] = df.apply(f2, axis=1)
        assert len(df[df['chr1']!=df['chr2']]) == 0
        df.rename(columns={'chr1': 'chr'}, inplace=True)
        df['x'] = round((df['x1'] + df['x2'])/2)
        df['y'] = round((df['y1'] + df['y2'])/2)
        df['cell_type'] = cell_type
        df = df[df.columns.intersection(['x', 'y', 'chr', 'cell_type'])]
        df = df.sort_values(by=['x'])
        new_dfs_dict[cell_type] = df

    return new_dfs_dict


def _concat_dfs(dfs_dict: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    """
    Concatenate a dataframes from a dictionary to one dataframe.
    Args:
        dictionary with partition ids as keys and pandas DataFrames as values.
    Returns:
        concatenated pandas DataFrame.
    """
    df = pd.concat(dfs_dict.values())
    df = df.sort_values(by=['x'])

    return df


def add_labels(dfs_dict: Dict[str, pd.DataFrame]) -> None:
    """
    Add labels to the dataframes depending on the cell type.
    Args:
        dictionary with partition ids as keys and pandas DataFrames as values.
    Returns:
        dictionary with partition ids as keys and changed pandas DataFrames as values.
    """
    df = _concat_dfs(dfs_dict)
    for name, cell_df in dfs_dict.items():
        cell_type = cell_df['cell_type'].unique()[0]
        f = lambda x: 1 if x['cell_type'] == cell_type else 0
        df['label'] = df.apply(f, axis=1)
        dfs_dict[name] = copy.deepcopy(df)
    
    return dfs_dict


def read_peaks(partitioned_input: Dict[str, Callable[[], Any]],
                cells2names: Dict[str, dict],
                dataset_name: str) -> Dict[str, pd.DataFrame]:
    """
    Load dataframes with DNase-seq/ChIP-seq peaks and modify the chromosome columns.
    Args:
        partitioned_input: dictionary with partition ids as keys and load functions as values.
        cells2names: dictionary, template: {'dataset_name': {'file_name': 'cell_type'}}
        dataset_name: name of the dataset to select from cells2names.
    Returns:
        dictionary with partition ids as keys and pandas DataFrames as values.
    """
    dfs_dict = _dict_partitions(partitioned_input)
    cells2names_dataset_dict = cells2names[dataset_name]
    keys_dict = {".".join(key.split(".")[:-1]): key for key in cells2names_dataset_dict}
    new_dfs_dict = dict()
    for name, df in dfs_dict.items():
        f = lambda x: x['chr'].split("chr")[-1]
        df["chr"] = df.apply(f, axis=1)
        df = df.sort_values(by=['start'])
        new_dfs_dict[cells2names_dataset_dict[keys_dict[name]]] = df

    return new_dfs_dict


def _count_peaks_single_df(main_df: pd.DataFrame, peaks_df: pd.DataFrame, experiment: str, r: int) -> pd.DataFrame:
    """
    Count the number of experiment peaks in both regions of each chromatin loop.
    Args:
        main_df: pandas DataFrame with chromatin loops.
        peaks_df: pandas DataFrame with experiment peaks.
    Returns:
        pandas DataFrame with chromatin loops and the numbers of experiment peaks in both regions of each loop
    """
    main_df['x_'+experiment+'_count'] = 0
    main_df['y_'+experiment+'_count'] = 0

    peaks_df = peaks_df.sort_values(by=['start'])

    # group by chromosome
    main_df_grouped = main_df.groupby('chr')
    peaks_df_grouped = peaks_df.groupby('chr')

    for main_chr, main_df_chr in main_df_grouped:
        for peaks_chr, peaks_df_chr in peaks_df_grouped:
            if main_chr == peaks_chr:
                    for idx, row in main_df_chr.iterrows():
                        for region in ['x', 'y']:
                            # TO CHANGE !!! - to overlapping peaks
                            # overlapping_peaks = peaks_df_chr[(peaks_df_chr['start'] >= row[region] - r)&
                            #                                     (peaks_df_chr['start'] <= row[region] + r)&
                            #                                     (peaks_df_chr['end'] >= row[region] - r)&
                            #                                     (peaks_df_chr['end'] <= row[region] + r)]
                            overlapping_peaks = peaks_df_chr[~(peaks_df_chr['start'] > row[region] + r)&
                                                             ~(peaks_df_chr['end'] < row[region] - r)]
                            main_df.loc[idx, f'{region}_{experiment}_count'] = len(overlapping_peaks)

    return main_df


def count_peaks(main_dfs_dict: Dict[str, Callable[[], Any]], peaks_dfs_dict: Dict[str, pd.DataFrame], experiment: str, r: int) -> pd.DataFrame:
    """
    Count the number of peaks in both regions of each chromatin loop,
    for each dataframe from the main_dfs_dict dictionary.
    Args:
        main_dfs_dict: dictionary with partition ids as keys and load functions of pandas DataFrames with chromatin loops as values.
        peaks_dfs_dict: dictionary with partition ids as keys and pandas DataFrames with experiment peaks as values.
        r: radius of the region around the loop center.
    Returns:
        dictionary with pandas DataFrames with chromatin loops and the numbers of peaks in both regions of each loop for each cell type.
    """
    print(f'Adding peaks counts for {experiment}...')
    # if main_dfs_dict values are not DataFrames, load them
    if not isinstance(list(main_dfs_dict.values())[0], pd.DataFrame):
        main_dfs_dict = _dict_partitions(main_dfs_dict)
    for peaks_name, peaks_df in tqdm(peaks_dfs_dict.items()):
        print(f'...for {peaks_name} cell...')
        for main_name, main_df in main_dfs_dict.items():
            if main_name == peaks_name:
                main_df = _count_peaks_single_df(main_df, peaks_df, experiment, r)
                main_dfs_dict[main_name] = main_df
    print('Done!')

    return main_dfs_dict


def read_bigWig(partitioned_input: Dict[str, Callable[[], Any]],
                cells2names: Dict[str, dict],
                dataset_name: str) -> Dict[str, pd.DataFrame]:
    """
    Load dataframes with DNase-seq/ChIP-seq bigWig data.
    Args:
        partitioned_input: dictionary with partition ids as keys and load functions as values.
        cells2names: dictionary, template: {'dataset_name': {'file_name': 'cell_type'}}
        dataset_name: name of the dataset to select from cells2names.
    Returns:
        dictionary with partition ids as keys and paths to bigWig files as values.
    """
    input_dict = _dict_partitions(partitioned_input)
    cells2names_dataset_dict = cells2names[dataset_name]
    keys_dict = {".".join(key.split(".")[:-1]): key for key in cells2names_dataset_dict}
    bigWig_data_dict = dict()
    for name, bigWig_path in input_dict.items():
        bigWig_data_dict[cells2names_dataset_dict[keys_dict[name]]] = bigWig_path

    return bigWig_data_dict


def _trimm_intervals(intervals: tuple, start: int, end: int):
    """
    Trim intervals to the range [start, end].
    Args:
        intervals: tuple with intervals.
        start: start of the range.
        end: end of the range.
    Returns:
        trimmed intervals.
    """
    trimmed_intervals = []
    for interval in intervals:
        if interval[0] < start:
            interval = (start, interval[1], interval[2])
        if interval[1] > end:
            interval = (interval[0], end, interval[2])
        trimmed_intervals.append(interval)

    return tuple(trimmed_intervals)


def _add_bigWig_data_single_df(main_df: pd.DataFrame, bigWig_path, experiment: str, r: int) -> pd.DataFrame:
    """
    Count the number of experiment peaks in both regions of each chromatin loop.
    Args:
        main_df: pandas DataFrame with chromatin loops.
        bigWig_obj: pyBigWig object with experiment peaks
    Returns:
        pandas DataFrame with chromatin loops and ...
    """
    bigWig_obj = pyBigWig.open(bigWig_path)
    regions = ['x', 'y']
    for region in regions:
        #main_df[f'{region}_{experiment}_values_vector'] = main_df.apply(lambda x: np.array(bigWig_obj.values('chr'+x['chr'], int(x[region]-r), int(x[region]+r)), dtype=np.float32), axis=1)
        main_df[f'{region}_{experiment}_mean'] = main_df.apply(lambda x: bigWig_obj.stats('chr'+x['chr'], int(x[region]-r), int(x[region]+r), type='mean')[0], axis=1)
        main_df[f'{region}_{experiment}_max'] = main_df.apply(lambda x: bigWig_obj.stats('chr'+x['chr'], int(x[region]-r), int(x[region]+r), type='max')[0], axis=1)
        main_df[f'{region}_{experiment}_min'] = main_df.apply(lambda x: bigWig_obj.stats('chr'+x['chr'], int(x[region]-r), int(x[region]+r), type='min')[0], axis=1)
    return main_df


def add_bigWig_data(main_dfs_dict: Dict[str, Callable[[], Any]],
                    bigWig_data_dict: dict,
                    experiment: str, r: int) -> pd.DataFrame:
    """
    Count the number of peaks in both regions of each chromatin loop,
    for each dataframe from the main_dfs_dict dictionary.
    Args:
        main_dfs_dict: dictionary with partition ids as keys and load functions of pandas DataFrames with chromatin loops as values.
        bigWig_data_dict: dictionary with cell types as keys and pyBigWig objects as values.
        r: radius of the region around the loop center.
    Returns:
        dictionary with pandas DataFrames with chromatin loops and ...
    """
    print(f'Adding {experiment} data...')
    # if main_dfs_dict values are not DataFrames, load them
    if not isinstance(list(main_dfs_dict.values())[0], pd.DataFrame):
        main_dfs_dict = _dict_partitions(main_dfs_dict)
    for bigWig_name, bigWig_path in tqdm(bigWig_data_dict.items()):
        print(f'...for {bigWig_name} cell...')
        for main_name, main_df in main_dfs_dict.items():
            if main_name == bigWig_name:
                main_df = _add_bigWig_data_single_df(main_df, bigWig_path, experiment, r)
                main_dfs_dict[main_name] = main_df
    print('Done!')

    return main_dfs_dict


def