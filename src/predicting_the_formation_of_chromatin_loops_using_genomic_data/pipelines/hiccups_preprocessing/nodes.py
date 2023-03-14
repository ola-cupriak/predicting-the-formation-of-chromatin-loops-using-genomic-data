from typing import Any, Callable, Dict
import pandas as pd
import copy



def _dict_partitions(partitioned_input: Dict[str, Callable[[], Any]]) -> Dict[str, pd.DataFrame]:
    """
    Load all partitions and save them in a dictionary.
    Args:
        partitioned_input: A dictionary with partition ids as keys and load functions as values.

    Returns:
        dictionary with partition ids as keys and pandas DataFrames as values.
    """
    result = dict()

    for partition_key, partition_load_func in sorted(partitioned_input.items()):
        partition_data = partition_load_func()  # load the actual partition data
        # concat with existing result
        result[partition_key] = partition_data

    return result


def preprocess_hiccups(partitioned_input: Dict[str, Callable[[], Any]]) -> Dict[str, pd.DataFrame]:
    """
    Modify the dataframes.
    Args:
        dictionary with partition ids as keys and pandas DataFrames as values.

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
        df['cell_type'] = name.split("_")[1]
        df = df[df.columns.intersection(['x', 'y', 'chr', 'cell_type'])]
        df = df.sort_values(by=['x'])
        dfs_dict[name] = df

    return dfs_dict


def _concat_dfs(dfs_dict: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    """
    Concatenate a dataframes from a dictionary.
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