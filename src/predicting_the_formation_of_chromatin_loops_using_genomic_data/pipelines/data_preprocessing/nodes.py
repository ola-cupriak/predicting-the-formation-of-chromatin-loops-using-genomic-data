import pyBigWig
from predicting_the_formation_of_chromatin_loops_using_genomic_data.utils import _dict_partitions
from typing import Any, Callable, Dict
import pandas as pd
import copy
from tqdm import tqdm
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import motifs
import pybedtools
import subprocess
import warnings

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
    df.loc[df['chr'] == 'X', 'chr'] = 100
    df.loc[df['chr'] == 'Y', 'chr'] = 200
    df['chr'] = df['chr'].astype(int)
    df = df.sort_values(by=['chr', region_name])
    df.loc[df['chr'] == 100, 'chr'] = 'X'
    df.loc[df['chr'] == 200, 'chr'] = 'Y'
    df['chr'] = df['chr'].astype(str)

    return df

    
def read_hic(partitioned_input: Dict[str, Callable[[], Any]], 
                cells2names: Dict[str, dict],
                dataset_name: str, r: int) -> Dict[str, pd.DataFrame]:
    """
    Load and modify the dataframes with chromatin loops anotations.
    Args:
        partitioned_input: dictionary with partition ids as keys and load functions as values.
        cells2names: dictionary, template: {'dataset_name': {'file_name': 'cell_type'}}
        dataset_name: name of the dataset to select from cells2names.
        r: radius of the region.
    Returns:
        dictionary with cell types as keys and pandas DataFrames as values.
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

        new_dfs_dict[cell_type] = df

    return new_dfs_dict


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
        dictionary with cell types as keys and pandas DataFrames as values.
    """
    dfs_dict = _dict_partitions(partitioned_input)
    cells2names_dataset_dict = cells2names[dataset_name]
    keys_dict = {".".join(key.split(".")[:-1]): key for key in cells2names_dataset_dict}
    new_dfs_dict = dict()
    for name, df in dfs_dict.items():
        f = lambda x: x['chr'].split("chr")[-1]
        df["chr"] = df.apply(f, axis=1)
        df = _sort_df(df, 'start')
        new_dfs_dict[cells2names_dataset_dict[keys_dict[name]]] = df

    return new_dfs_dict


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
        dictionary with cell types as keys and paths to bigWig files as values.
    """
    input_dict = _dict_partitions(partitioned_input)
    cells2names_dataset_dict = cells2names[dataset_name]
    keys_dict = {".".join(key.split(".")[:-1]): key for key in cells2names_dataset_dict}
    bigWig_data_dict = dict()
    for name, bigWig_path in input_dict.items():
        bigWig_data_dict[cells2names_dataset_dict[keys_dict[name]]] = bigWig_path

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


def add_labels(dfs_dict: Dict[str, pd.DataFrame]) -> None:
    """
    Add labels to the dataframes depending on the cell type.
    Args:
        dfs_dict: dictionary with cell types as keys and pandas DataFrames as values.
    Returns:
        dictionary with cell types as keys and changed pandas DataFrames as values.
    """
    df = _concat_dfs(dfs_dict)
    df = _sort_df(df, 'x')
    for name, cell_df in dfs_dict.items():
        cell_type = cell_df['cell_type'].unique()[0]
        f = lambda x: 1 if x['cell_type'] == cell_type else 0
        df['label'] = df.apply(f, axis=1)
        dfs_dict[name] = copy.deepcopy(df)
    
    return dfs_dict


def get_overlapping_regions(df1: pd.DataFrame, df2: pd.DataFrame, count: bool = False) -> pd.DataFrame:
    """
    Get overlapping regions from two bed files.
    Args:
        df1: path to bed file 1
        df2: path to bed file 2
        count: if True, count the number of overlapping regions.
    Returns:
        pd.DataFrame with overlapping regions or number of overlapping regions.
    """
    bed1 = pybedtools.BedTool.from_dataframe(df1)
    bed2 = pybedtools.BedTool.from_dataframe(df2)
    if count:
        intersection = bed1.intersect(bed2, c=count)
        intersection = pd.read_table(intersection.fn, names=['chr', 'start', 'end', 'second_reg', 'cell_type', 'count'])
    else:
        intersection = bed1.intersect(bed2)
        intersection = pd.read_table(intersection.fn, names=['chr', 'start', 'end'])
    
    

    return intersection


def _count_peaks_single_df(main_df: pd.DataFrame, peaks_df: pd.DataFrame, experiment: str, r: int) -> pd.DataFrame:
    """
    Count the number of experiment peaks in both regions of each chromatin loop.
    Args:
        main_df: pandas DataFrame with chromatin loops.
        peaks_df: pandas DataFrame with experiment peaks.
        anchors_df: pandas DataFrame with anchors coordinates and sequences.
        experiment: name of the experiment.
        r: radius of the anchors.
        motifs_path: path to the file with motifs.
    Returns:
        pandas DataFrame with chromatin loops and the numbers of experiment peaks in both regions of each loop
        and columns with counts of each motif in each region if motifs_path is not None.
    """
    for region in ['x', 'y']:
        second_reg = 'y' if region == 'x' else 'x'
        peak_counts = main_df[['chr', f'{region}_start', f'{region}_end', f'{second_reg}', 'cell_type']]
        peak_counts = get_overlapping_regions(peak_counts, peaks_df, count=True)
        peak_counts.rename(columns={'count': f'{region}_{experiment}_counts',
                                    'start': f'{region}_start',
                                    'end': f'{region}_end',
                                    'second_reg': f'{second_reg}'}, inplace=True)

        assert len(peak_counts) == len(main_df), 'Length of the main_df and peak_counts are not equal'

        main_len_before = len(main_df)
        main_df = main_df.merge(peak_counts)

        assert len(main_df) == main_len_before, 'Length of the main_df changed after merging'

    return main_df



def count_peaks(main_dfs_dict: Dict[str, Callable[[], Any]], peaks_dfs_dict: Dict[str, pd.DataFrame], 
                experiment: str, r: int) -> pd.DataFrame:
    """
    Count the number of peaks in both regions of each chromatin loop,
    for each dataframe from the main_dfs_dict dictionary.
    Args:
        main_dfs_dict: dictionary with cell types as keys and load functions of pandas DataFrames with chromatin loops as values.
        peaks_dfs_dict: dictionary with cell types as keys and pandas DataFrames with experiment peaks as values.
        anchors_df: pandas DataFrame with anchors coordinates and sequences.
        experiment: name of the experiment.
        r: radius of the region around the loop center.
        motifs_path: path to the file with motifs.
    Returns:
        dictionary with cell types as keys and pandas DataFrames with chromatin loops and the numbers of peaks 
        (and motifs counts if motifs_path is not None) in both regions of each loop as values.
    """
    print(f'Adding peaks counts for {experiment}...')
    # if main_dfs_dict values are not DataFrames, load them
    if not isinstance(list(main_dfs_dict.values())[0], pd.DataFrame):
        main_dfs_dict = _dict_partitions(main_dfs_dict)
    for peaks_name, peaks_df in peaks_dfs_dict.items():
        print(f'...for {peaks_name} cell...')
        for main_name, main_df in main_dfs_dict.items():
            if main_name == peaks_name:
                main_df = _count_peaks_single_df(main_df, peaks_df, experiment, r)
                main_dfs_dict[main_name] = main_df
    print('Done!')

    return main_dfs_dict


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
        main_df[f'{region}_{experiment}_mean'] = main_df.apply(lambda x: bigWig_obj.stats('chr'+x['chr'], x[f'{region}_start'], x[f'{region}_end'], type='mean')[0], axis=1)
        main_df[f'{region}_{experiment}_max'] = main_df.apply(lambda x: bigWig_obj.stats('chr'+x['chr'], x[f'{region}_start'], x[f'{region}_end'], type='max')[0], axis=1)
        main_df[f'{region}_{experiment}_min'] = main_df.apply(lambda x: bigWig_obj.stats('chr'+x['chr'], x[f'{region}_start'], x[f'{region}_end'], type='min')[0], axis=1)
    return main_df


def add_bigWig_data(main_dfs_dict: Dict[str, Callable[[], Any]],
                    bigWig_data_dict: dict,
                    experiment: str, r: int) -> pd.DataFrame:
    """
    Count the number of peaks in both regions of each chromatin loop,
    for each dataframe from the main_dfs_dict dictionary.
    Args:
        main_dfs_dict: dictionary with cell types as keys and load functions of pandas DataFrames with chromatin loops as values.
        bigWig_data_dict: dictionary with cell types as keys and pyBigWig objects as values.
        r: radius of the region around the loop center.
    Returns:
        dictionary with pandas DataFrames with chromatin loops and ...
    """
    print(f'Adding {experiment} data...')
    # if main_dfs_dict values are not DataFrames, load them
    if not isinstance(list(main_dfs_dict.values())[0], pd.DataFrame):
        main_dfs_dict = _dict_partitions(main_dfs_dict)
    for bigWig_name, bigWig_path in bigWig_data_dict.items():
        print(f'...for {bigWig_name} cell...')
        for main_name, main_df in main_dfs_dict.items():
            if main_name == bigWig_name:
                main_df = _add_bigWig_data_single_df(main_df, bigWig_path, experiment, r)
                main_dfs_dict[main_name] = main_df
    print('Done!')

    return main_dfs_dict


def all_anchors2one_df(dfs_dict: Dict[str, pd.DataFrame], r: int) -> pd.DataFrame:
    """
    Combines two columns with regions into one colun and removes duplicate rows.
    Args:
        dfs_dict: dictionary with cell types as keys and pandas DataFrames as values.
        r: radius of the region.
    Returns:
        pandas DataFrame with one column with regions.
    """
    df_with_2_regions = _concat_dfs(dfs_dict)
    df_part1 = df_with_2_regions[['chr', 'x_start', 'x_end']]
    df_part2 = df_with_2_regions[['chr', 'y_start', 'y_end']]
    # rename columns
    df_part1.columns = ['chr', 'start', 'end']
    df_part2.columns = ['chr', 'start', 'end']
    # combine two columns into one
    anchors_df = pd.concat([df_part1, df_part2], axis=0)
    # remove duplicate rows
    anchors_df = anchors_df.drop_duplicates()
    # sort by chr and region
    anchors_df = _sort_df(anchors_df, 'start')
    # reset index
    anchors_df = anchors_df.reset_index(drop=True)
    anchors_df = anchors_df[['chr', 'start', 'end']]
    
    return anchors_df


def all_peaks2one_df(peaks_dict: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    df = _concat_dfs(peaks_dict)
    # sort by chr and region
    df = _sort_df(df, 'start')
    # remove duplicate rows
    df = df.drop_duplicates()
    # reset index
    df = df.reset_index(drop=True)
    df = df[['chr', 'start', 'end']]

    return df
        

def getfasta_bedfile(df: pd.DataFrame, path_simp_genome: str) -> str:
    """
    Cut sequences from chromosomes using bedtools for coordinates from bed file.
    Args:
        df: pandas DataFrame with coordinates.
        path_simp_genome: path to fasta file with chromosomes.
    Returns:
        string with fasta sequences.
    """
    fasta = pybedtools.BedTool(path_simp_genome)
    df_to_search = df
    df_to_search['chr'] = 'chr' + df_to_search['chr']
    bed = pybedtools.BedTool.from_dataframe(df_to_search)
    fasta_bed = bed.sequence(fi=fasta)

    records = []
    for record in SeqIO.parse(fasta_bed.seqfn, 'fasta'):
        records.append(record)

    return records


def anchors2df(anchors_fasta: list) -> pd.DataFrame:
    """Convert anchors fasta file to pandas DataFrame.
    Args:
        anchors_fasta (str): string with fasta sequences for anchors.
    Returns:
        pandas DataFrame with anchors coordinates and sequences.
    """
    dict2df = {'chr': [], 'start': [], 'end': [], 'seq': []}
    for record in anchors_fasta:
        dict2df['chr'].append(record.id.split(':')[0].split('chr')[-1])
        dict2df['start'].append(int(record.id.split(':')[1].split('-')[0]))
        dict2df['end'].append(int(record.id.split(':')[1].split('-')[1]))
        dict2df['seq'].append(str(record.seq))

    anchors_df = pd.DataFrame(dict2df)
    anchors_df = anchors_df.sort_values(by=['chr', 'start'])

    return anchors_df