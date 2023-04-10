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
from io import StringIO
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
        # set dtypes
        df = df.astype({'x': 'int32', 'y': 'int32', 'x_start': 'int32', 'x_end': 'int32', 
                        'y_start': 'int32', 'y_end': 'int32', 'chr': 'string', 'cell_type': 'string'})

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
        cell_type = cells2names_dataset_dict[keys_dict[name]]
        f = lambda x: x['chr'].split("chr")[-1]
        df["chr"] = df.apply(f, axis=1)
        df = _sort_df(df, 'start')
        df['cell_type'] = cell_type
        df['cell_type'] = df['cell_type'].astype('string')
        new_dfs_dict[cell_type] = df

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
        df['label'] = df['label'].astype('int16')
        dfs_dict[name] = copy.deepcopy(df)
    
    return dfs_dict


def get_overlapping_regions(df1: pd.DataFrame, df2: pd.DataFrame, count: bool = False, wa: bool = False, wb: bool = False) -> pd.DataFrame:
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
        if wa and wb:
            intersection = bed1.intersect(bed2, wa=wa, wb=wb)
            intersection = pd.read_table(intersection.fn, names=['anchor_chr', 'anchor_start', 'anchor_end', 'peak_chr', 'peak_start', 'peak_end', 'cell_type'])
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

    to_change_dtype = ['x_DNase_seq_peaks_counts', 'y_DNase_seq_peaks_counts']
    main_df[to_change_dtype] = main_df[to_change_dtype].astype('int16')

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
        # change dtype to save memory
        to_change_dtype = [f'{region}_{experiment}_mean', f'{region}_{experiment}_max', f'{region}_{experiment}_min']
        main_df[to_change_dtype] = main_df[to_change_dtype].astype('float32')
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
    df = df.reset_index(drop=True)

    return df


def _prepare_peaks_df(peaks_df: pd.DataFrame):
    
    peaks_df = peaks_df[['chr', 'start', 'end']]
    # remove duplicate rows
    peaks_df = peaks_df.drop_duplicates()
    # reset index
    peaks_df = peaks_df.reset_index(drop=True)

    return peaks_df


def get_overlaps_with_names(anchors_df: pd.DataFrame, peaks_df: pd.DataFrame) -> pd.DataFrame:
    """
    """
    #peaks_df = _prepare_peaks_df(peaks_df)

    intersection_wa_wb = get_overlapping_regions(anchors_df, peaks_df, wa=True, wb=True)
    intersection = get_overlapping_regions(anchors_df, peaks_df)

    joined_intersection = pd.merge(intersection_wa_wb, intersection, left_index=True, right_index=True)
    joined_intersection['name'] = joined_intersection.apply(lambda x: f"{x['chr']}:{x['anchor_start']}:{x['anchor_end']}:{x['cell_type']}", axis=1)

    return joined_intersection[['chr', 'start', 'end', 'name']]


def getfasta_bedfile(df: pd.DataFrame, path_simp_genome: str, path_to_save: str) -> str:
    """
    Cut sequences from chromosomes using bedtools for coordinates from bed file.
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


def _motify_output(lines: list):
    motif_id2name = {}
    for i in tqdm(range(len(lines))):
        if i == 0:
            lines[i] = lines[i].replace('sequence_name', 'chr\tstart\tend\tcell_type')
        else:
            elements = lines[i].split('\t')
            motif_id, motif_alt_id, name, _, _, strand, _, _, _, _ = elements

            motif_id2name[motif_id] = motif_alt_id

            name = name.replace(':', '\t')
            if strand == '+':
                motif_id += '_f'
            else:
                motif_id += '_r'
            
            elements[0] = motif_id
            elements[2] = name

            lines[i] = '\t'.join(elements)
    
    with open('data/02_intermediate/motif_id2name.txt', 'w') as file:
        for key, value in motif_id2name.items():
            file.write(f'{key}\t{value}\n')
    
    return lines



def find_motifs(path_motifs: str, path_fasta: list) -> pd.DataFrame:
    """
    Find motifs in fasta sequences.
    Args:
        path_motifs: path to file with motifs in jaspar format.
        path_fasta: path to fasta file with sequences.
    Returns:
        pandas DataFrame with counts of motifs found in anchors.
    """
    # change jaspar format to meme format
    path_for_meme = path_motifs.replace('.txt', '.meme')
    subprocess.run(f'jaspar2meme -bundle {path_motifs} > {path_for_meme}', shell=True)
    # find motifs
    start = time.time()
    print('Finding motifs...')
    subprocess.run(f'fimo --text {path_for_meme} {path_fasta} > data/temp/temp.csv', shell=True)
    print(f'Finding motifs took {time.time() - start} seconds')
    # read output
    with open('data/temp/temp.csv', 'r') as f:
        output = f.readlines()
    # save modified output to temp file
    output = _motify_output(output)
    with open('data/temp/temp.csv', 'w') as f:
        f.writelines(output)
    output = None
    # read temp file as pandas DataFrame
    dtypes = {'chr': "string", 'start': "int32", 'end': "int32", 'motif_id': "string", 'cell_type': "string"}
    df = pd.read_csv('data/temp/temp.csv', sep='\t', dtype=dtypes, usecols=list(dtypes.keys()))
    subprocess.run('rm data/temp/temp.csv', shell=True)
    df = df[['chr', 'start', 'end', 'motif_id', 'cell_type']]
    len_before = len(df)
    # count motif occurences
    df = pd.DataFrame(df.groupby(['chr', 'start', 'end', 'cell_type', 'motif_id'], observed=True).size(), columns=['count'])
    df['count'] = df['count'].astype('int16')
    df = df.unstack('motif_id', fill_value=0)
    df.columns = ['_'.join(x) for x in df.columns if x[0]=='count']
    df.reset_index(inplace=True)
    df.rename(columns={name: name.replace('count_', '') for name in df.columns}, inplace=True)
    
    assert sum(df.iloc[:, 4:].sum(axis=1)) == len_before, 'Something went wrong with counting motifs'
    
    return df


def _reverse_names(colnames: list):
    
    to_reverse_f = [name for name in colnames if '_f' in name]
    to_reverse_r = [name for name in colnames if '_r' in name]
    to_reverse_f = {name: name.replace('_f', '_r') for name in to_reverse_f}
    to_reverse_r = {name: name.replace('_r', '_f') for name in to_reverse_r}

    to_reverse = {**to_reverse_f, **to_reverse_r}

    return to_reverse


def _count_motifs_single_df(main_df: pd.DataFrame, motifs_df: pd.DataFrame) -> pd.DataFrame:
    
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