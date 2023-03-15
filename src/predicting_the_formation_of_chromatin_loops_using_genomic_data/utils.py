from typing import Any, Callable, Dict
import pandas as pd



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