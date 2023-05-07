import copy
import numpy as np
import pyBigWig
from predicting_the_formation_of_chromatin_loops_using_genomic_data.utils import _dict_partitions
from typing import Any, Callable, Dict
import pandas as pd
from tqdm import tqdm
from Bio import SeqIO
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
        cell_type = cells2names_dataset_dict[keys_dict[name]]
        if cells_to_use and cell_type not in cells_to_use:
            continue
        f = lambda x: x['chr'].split("chr")[-1]
        df["chr"] = df.apply(f, axis=1)
        df = _sort_df(df, 'start')
        df['cell_type'] = cell_type
        df['cell_type'] = df['cell_type'].astype('string')
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


def add_labels(dfs_dict: Dict[str, pd.DataFrame]) -> None:
    """
    Add labels to the dataframes depending on the cell type.
    Args:
        dfs_dict: dictionary with cell types as keys and pandas DataFrames as values.
    Returns:
        dictionary:
            keys: cell types
            values: pandas DataFrames with labels added.
    """
    df = _concat_dfs(dfs_dict)
    df = _sort_df(df, 'x')
    for name, cell_df in dfs_dict.items():
        cell_type = cell_df['cell_type'].unique()[0]
        f = lambda x: 1 if x['cell_type'] == cell_type else 0
        df['label'] = df.apply(f, axis=1)
        df['label'] = df['label'].astype('int16')
        dfs_dict[name] = copy.deepcopy(df)
    
    return dfs_dict


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
    for region in regions:
        main_df[f'{region}_{experiment}_weighted_mean'] = main_df.apply(lambda x: calculate_weighted_mean(bigWig_obj.values('chr'+x['chr'], x[f'{region}_start'], x[f'{region}_end'])), axis=1)
        main_df[f'{region}_{experiment}_mean'] = main_df.apply(lambda x: bigWig_obj.stats('chr'+x['chr'], x[f'{region}_start'], x[f'{region}_end'], type='mean')[0], axis=1)
        main_df[f'{region}_{experiment}_max'] = main_df.apply(lambda x: bigWig_obj.stats('chr'+x['chr'], x[f'{region}_start'], x[f'{region}_end'], type='max')[0], axis=1)
        main_df[f'{region}_{experiment}_min'] = main_df.apply(lambda x: bigWig_obj.stats('chr'+x['chr'], x[f'{region}_start'], x[f'{region}_end'], type='min')[0], axis=1)
        # change dtype to save memory
        to_change_dtype = [f'{region}_{experiment}_weighted_mean', f'{region}_{experiment}_mean', f'{region}_{experiment}_max', f'{region}_{experiment}_min']
        main_df[to_change_dtype] = main_df[to_change_dtype].astype('float32')
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
    Concatenates all anchors DataFrames into one DataFrame and 
    combines the columns describing x regions and the columns describing y regions 
    into one set of columns describing all regions.
    (columns: x_chr, x_start, x_end, y_chr, y_start, y_end -> columns: chr, start, end)
    Args:
        dfs_dict: dictionary with cell types as keys and pandas DataFrames as values.
    Returns:
        pandas DataFrame with one set of columns describing all regions.
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


def _modify_output(lines: list):
    """
    Modifies each line of output of the FIMO tool.
    Saves motif IDs and corresponding names to a file.
    Args:
        lines: list of lines from the output file.
    Returns:
        list of modified lines.
    """
    motif_id2name = {}
    for i in tqdm(range(len(lines))):
        # Change header to tsv format
        if i == 0:
            lines[i] = lines[i].replace('sequence_name', 'chr\tstart\tend\tcell_type')
        else:
            # Get the information you need 
            elements = lines[i].split('\t')
            motif_id, motif_alt_id, name, _, _, strand, _, _, _, _ = elements
            # Save motif_id and motif_alt_id to the dictionary
            motif_id2name[motif_id] = motif_alt_id
            # Change name format to tsv format
            name = name.replace(':', '\t')
            # Add strand information to the motif_id
            if strand == '+':
                motif_id += '_f'
            else:
                motif_id += '_r'
            # Update the list of lines
            elements[0] = motif_id
            elements[2] = name
            lines[i] = '\t'.join(elements)
    # Save motif_id2name to the file
    with open('data/02_intermediate/motif_id2name.txt', 'w') as file:
        for key, value in motif_id2name.items():
            file.write(f'{key}\t{value}\n')
    
    return lines


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
    with open('data/temp/temp.csv', 'r') as f:
        output = f.readlines()
    # Save modified output to temporary file
    output = _modify_output(output)
    with open('data/temp/temp.csv', 'w') as f:
        f.writelines(output)
    output = None
    # Read temporary file as pandas DataFrame
    dtypes = {'chr': "string", 'start': "int32", 'end': "int32", 'motif_id': "string", 'cell_type': "string"}
    df = pd.read_csv('data/temp/temp.csv', sep='\t', dtype=dtypes, usecols=list(dtypes.keys()))
    subprocess.run('rm data/temp/temp.csv', shell=True)
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
    assert len(main_df) == len(main_df['id'].unique()), 'Something went wrong with adding id column'
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
    x_overlap = _get_overlapping_regions(df_1x, df_0x, names=['chr1', 'satrt1', 'end1', 'chr2', 'start2', 'end2', 'id'], wa=True, wb=True)
    y_overlap = _get_overlapping_regions(df_1y, df_0y, names=['chr1', 'satrt1', 'end1', 'chr2', 'start2', 'end2', 'id'], wa=True, wb=True)
    # Find overlapping regions for x and y anchors together
    x_and_y = (set(x_overlap['id']) & set(y_overlap['id']))
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