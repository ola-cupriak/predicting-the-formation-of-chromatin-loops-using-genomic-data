from predicting_the_formation_of_chromatin_loops_using_genomic_data.utils import _dict_partitions
from typing import Any, Callable, Dict
import pandas as pd
import copy



def read_hiccups(partitioned_input: Dict[str, Callable[[], Any]]) -> Dict[str, pd.DataFrame]:
    """
    Load and modify the dataframes with chromatin loops anotations.
    Args:
        dictionary with partition ids as keys and load functions as values.

    Returns:
        dictionary with partition ids as keys and pandas DataFrames as values.
    """
    dfs_dict = _dict_partitions(partitioned_input)
    for name, df in dfs_dict.items():
        f1 = lambda x: x['chr1'].split("chr")[-1]
        f2 = lambda x: x['chr2'].split("chr")[-1]
        df["chr1"] = df.apply(f1, axis=1)
        df["chr2"] = df.apply(f2, axis=1)
        assert len(df[df['chr1']!=df['chr2']]) == 0
        df.rename(columns={'chr1': 'chr'}, inplace=True)
        df['x'] = round((df['x1'] + df['x2'])/2)
        df['y'] = round((df['y1'] + df['y2'])/2)
        df['cell_type'] = name.split("_")[0]
        df = df[df.columns.intersection(['x', 'y', 'chr', 'cell_type'])]
        df = df.sort_values(by=['x'])
        dfs_dict[name] = df

    return dfs_dict


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


def read_signals(partitioned_input: Dict[str, Callable[[], Any]]) -> Dict[str, pd.DataFrame]:
    """
    Load dataframes with DNase-seq/ChIP-seq signals and modify the chromosome columns.
    Args:
        dictionary with partition ids as keys and load functions as values.
    Returns:
        dictionary with partition ids as keys and pandas DataFrames as values.
    """
    dfs_dict = _dict_partitions(partitioned_input)
    for name, df in dfs_dict.items():
        f = lambda x: x['chr'].split("chr")[-1]
        df["chr"] = df.apply(f, axis=1)
        df = df.sort_values(by=['start'])
        dfs_dict[name] = df

    return dfs_dict


def _count_signals(main_df: pd.DataFrame, signals_df: pd.DataFrame, experiment: str, r: int) -> pd.DataFrame:
    """
    Count the number of experiment signals in both regions of each chromatin loop.
    Args:
        main_df: pandas DataFrame with chromatin loops.
        signals_df: pandas DataFrame with experiment signals.
    Returns:
        pandas DataFrame with chromatin loops and the numbers of experiment signals in both regions of each loop
    """
    main_df[experiment+'_x_signals_count'] = 0
    main_df[experiment+'_y_signals_count'] = 0

    signals_df = signals_df.sort_values(by=['start'])

    # group by chromosome
    main_df_grouped = main_df.groupby('chr')
    signals_df_grouped = signals_df.groupby('chr')

    for main_chr, main_df_chr in main_df_grouped:
        for signals_chr, signals_df_chr in signals_df_grouped:
            if main_chr == signals_chr:
                    for idx, row in main_df_chr.iterrows():
                        for region in ['x', 'y']:
                            overlapping_signals = signals_df_chr[(signals_df_chr['start'] >= row[region] - r)&
                                                                (signals_df_chr['start'] <= row[region] + r)&
                                                                (signals_df_chr['end'] >= row[region] - r)&
                                                                (signals_df_chr['end'] <= row[region] + r)]
                            main_df.loc[idx, f'{experiment}_{region}_signals_count'] = len(overlapping_signals)

    return main_df


def count_signals(main_dfs_dict: Dict[str, Callable[[], Any]], signals_dfs_dict: Dict[str, pd.DataFrame], experiment: str, r: int) -> pd.DataFrame:
    """
    Count the number of signals in both regions of each chromatin loop,
    for each dataframe from the main_dfs_dict dictionary.
    Args:
        main_dfs_dict: dictionary with partition ids as keys and load functions of pandas DataFrames with chromatin loops as values.
        signals_dfs_dict: dictionary with partition ids as keys and pandas DataFrames with experiment signals as values.
        r: radius of the region around the loop center.
    Returns:
        pandas DataFrame with chromatin loops and the numbers of signals in both regions of each loop
    """
    main_dfs_dict = _dict_partitions(main_dfs_dict)
    for signals_name, signals_df in signals_dfs_dict.items():
        for main_name, main_df in main_dfs_dict.items():
            if main_name.split('_')[0] == signals_name.split('_')[0]:
                main_df = _count_signals(main_df, signals_df, experiment, r)
                main_dfs_dict[main_name] = main_df

    return main_dfs_dict

