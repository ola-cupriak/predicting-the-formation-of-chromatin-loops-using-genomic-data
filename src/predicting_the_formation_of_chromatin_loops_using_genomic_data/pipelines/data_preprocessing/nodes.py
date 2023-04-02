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
import subprocess
import warnings

warnings.simplefilter(action="ignore")


def simplify_genome_file(path: str) -> list:
    """Get chromosomes dict from fasta file.
    Args:
        path (str): path to fasta file with genome.
    Returns:
        list of SeqRecord objects.
    """
    genome = SeqIO.parse(open(path), 'fasta')
    chromosomes = []
    for sequence in genome:
        desc = sequence.description
        seq_id = sequence.id
        if seq_id.startswith('NC') and 'chromosome' in desc:
            chrom = 'chr' + desc.split('chromosome ')[1].split(',')[0]
            chromosomes.append(SeqRecord(sequence.seq, id=chrom, description=''))
    
    return chromosomes

    
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
        df['x'] = round((df['x1'] + df['x2'])/2)
        df['y'] = round((df['y1'] + df['y2'])/2)
        df['cell_type'] = cell_type
        df = df[df.columns.intersection(['x', 'y', 'chr', 'cell_type'])]
        df = df.sort_values(by=['x'])
        if cell_type == 'HMEC':
            new_dfs_dict[cell_type] = df

    return new_dfs_dict


def _concat_dfs(dfs_dict: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    """
    Concatenate a dataframes from a dictionary to one dataframe.
    Args:
        dictionary with cell types as keys and pandas DataFrames as values.
    Returns:
        concatenated pandas DataFrame.
    """
    df = pd.concat(list(dfs_dict.values()))
    df = df.sort_values(by=['x'])

    return df


def gather_all_anchors_into_df(dfs_dict: Dict[str, pd.DataFrame], r: int) -> pd.DataFrame:
    """
    Combines two columns with regions into one colun and removes duplicate rows.
    Args:
        dfs_dict: dictionary with cell types as keys and pandas DataFrames as values.
    Returns:
        pandas DataFrame with one column with regions.
    """
    df2cols = _concat_dfs(dfs_dict)
    df_part1 = df2cols[['chr', 'x']]
    df_part2 = df2cols[['chr', 'y']]
    # rename columns
    df_part1.columns = ['chr', 'anchor']
    df_part2.columns = ['chr', 'anchor']
    # combine two columns into one
    anchors_df = pd.concat([df_part1, df_part2], axis=0)
    # remove duplicate rows
    anchors_df = anchors_df.drop_duplicates()
    # add start and end columns
    anchors_df['start'] = (anchors_df['anchor'] - r).astype(int)
    anchors_df['end'] = (anchors_df['anchor'] + r).astype(int)
    # sort by chr and region
    anchors_df = anchors_df.sort_values(by=['chr', 'start'])
    anchors_df['chr'] = 'chr' + anchors_df['chr'].astype(str)
    # reset index
    anchors_df = anchors_df.reset_index(drop=True)
    anchors_df = anchors_df[['chr', 'start', 'end']]

    return anchors_df


def getfasta_anchors(path_anchors: str, path_simp_genome: str) -> str:
    """
    Cut sequences for anchors from chromosomes using bedtools.
    Args:
        path_anchors: path to bed file with anchors.
        path_simp_genome: path to fasta file with chromosomes.
    Returns:
        string with fasta sequences for anchors.
    """
    proc = subprocess.Popen(f'bedtools getfasta -fi {path_simp_genome} -bed {path_anchors}', shell=True, stdout=subprocess.PIPE)
    output = proc.stdout.read().decode("utf-8").split('>')[1:]

    records = []

    for record in output:
        record_id = record.split('\n')[0]
        record_seq = ''.join(record.split('\n')[1:])
        records.append(SeqRecord(Seq(record_seq), id=record_id, description=''))

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


def add_labels(dfs_dict: Dict[str, pd.DataFrame]) -> None:
    """
    Add labels to the dataframes depending on the cell type.
    Args:
        dfs_dict: dictionary with cell types as keys and pandas DataFrames as values.
    Returns:
        dictionary with cell types as keys and changed pandas DataFrames as values.
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
        dictionary with cell types as keys and pandas DataFrames as values.
    """
    dfs_dict = _dict_partitions(partitioned_input)
    cells2names_dataset_dict = cells2names[dataset_name]
    keys_dict = {".".join(key.split(".")[:-1]): key for key in cells2names_dataset_dict}
    new_dfs_dict = dict()
    for name, df in dfs_dict.items():
        f = lambda x: x['chr'].split("chr")[-1]
        df["chr"] = df.apply(f, axis=1)
        df = df.sort_values(by=['start'])
        if cells2names_dataset_dict[keys_dict[name]] == 'HMEC':
            new_dfs_dict[cells2names_dataset_dict[keys_dict[name]]] = df

    return new_dfs_dict



def _count_peaks_single_df(main_df: pd.DataFrame, peaks_df: pd.DataFrame, 
                           anchors_df: pd.DataFrame, experiment: str, r: int, 
                           motifs_path: str=None,) -> pd.DataFrame:
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
    # add columns for peak counts
    main_df['x_'+experiment+'_count'] = 0
    main_df['y_'+experiment+'_count'] = 0
    
    if motifs_path:
        motifs_dict = {}
        motifs_list = []
        for m in motifs.parse(open(motifs_path), "jaspar"):
            # create unique names for each motif
            motif_name = m.name.upper()
            if motif_name not in motifs_dict.keys():
                motifs_dict[motif_name] = 0
            else:
                motifs_dict[motif_name] += 1
                motif_name = f'{motif_name}_{motifs_dict[motif_name]}'
            motifs_list.append((motif_name, m.pssm))
            assert f'x_{motif_name}_f' not in main_df.columns, f'column {motif_name} already exists in main_df'

            # add columns for each motif counts
            main_df[f'x_motif_{motif_name}_f'] = 0
            main_df[f'x_motif_{motif_name}_b'] = 0
            main_df[f'y_motif_{motif_name}_f'] = 0
            main_df[f'y_motif_{motif_name}_b'] = 0

    # group by chromosome
    main_df_grouped = main_df.groupby('chr')
    peaks_df_grouped = peaks_df.groupby('chr')

    for main_chr, main_df_chr in tqdm(main_df_grouped):
        for peaks_chr, peaks_df_chr in peaks_df_grouped:
            # if chromosome is the same
            if main_chr == peaks_chr:
                    for idx, row in main_df_chr.iterrows():
                        for region in ['x', 'y']:
                            anchor_start = row[region] - r
                            anchor_end = row[region] + r
                            # find all peaks overlapping with anchor and count them
                            overlapping_peaks = peaks_df_chr[~(peaks_df_chr['start'] > anchor_end)&
                                                             ~(peaks_df_chr['end'] < anchor_start)]
                            main_df.loc[idx, f'{region}_{experiment}_count'] = len(overlapping_peaks)
                            
                            # find motif counts in anchor
                            if motifs_path:
                                # find anchor sequence
                                anchor_record = anchors_df[(anchors_df['start'] == anchor_start)&(anchors_df['chr'] == main_chr)]
                                assert len(anchor_record) == 1, f'not found or more than 1 anchor with coordinates found: {main_chr}:{anchor_start}-{anchor_end}'
                                anchor_seq = anchor_record.iloc[0]['seq']
                                # iterate over motifs
                                for motif_name, motif_pssm in motifs_list:
                                    count_peak_f = 0
                                    count_peak_b = 0
                                    # iterate over peaks of open chromatin
                                    for _, peak in overlapping_peaks.iterrows():
                                        peak_start = max(peak['start'], anchor_start)
                                        peak_end = min(peak['end'], anchor_end)
                                        seq_to_search = anchor_seq[int(peak_start-anchor_start):int(peak_end-anchor_start+1)]
                                        # check if motif is not longer than region of open chromatin to search
                                        if not len(motif_pssm[0]) > len(seq_to_search):
                                            for position, _ in motif_pssm.search(seq_to_search):
                                                if position >= 0:
                                                    if region == 'x':
                                                        count_peak_f += 1
                                                    else:
                                                        count_peak_b += 1
                                                else:
                                                    if region == 'x':
                                                        count_peak_b += 1
                                                    else:
                                                        count_peak_f += 1

                                    main_df.loc[idx, f'{region}_motif_{motif_name}_f'] = count_peak_f
                                    main_df.loc[idx, f'{region}_motif_{motif_name}_b'] = count_peak_b

    return main_df


def count_peaks(main_dfs_dict: Dict[str, Callable[[], Any]], peaks_dfs_dict: Dict[str, pd.DataFrame], 
                anchors_df: pd.DataFrame, experiment: str, r: int, motifs_path: str=None) -> pd.DataFrame:
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
    if motifs_path:
        print(f'Adding peaks counts for {experiment} and motifs counts...')
    else:
        print(f'Adding peaks counts for {experiment}...')
    # if main_dfs_dict values are not DataFrames, load them
    if not isinstance(list(main_dfs_dict.values())[0], pd.DataFrame):
        main_dfs_dict = _dict_partitions(main_dfs_dict)
    for peaks_name, peaks_df in peaks_dfs_dict.items():
        print(f'...for {peaks_name} cell...')
        for main_name, main_df in main_dfs_dict.items():
            if main_name == peaks_name:
                main_df = _count_peaks_single_df(main_df, peaks_df, anchors_df, experiment, r, motifs_path)
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
        dictionary with cell types as keys and paths to bigWig files as values.
    """
    input_dict = _dict_partitions(partitioned_input)
    cells2names_dataset_dict = cells2names[dataset_name]
    keys_dict = {".".join(key.split(".")[:-1]): key for key in cells2names_dataset_dict}
    bigWig_data_dict = dict()
    for name, bigWig_path in input_dict.items():
        bigWig_data_dict[cells2names_dataset_dict[keys_dict[name]]] = bigWig_path

    return bigWig_data_dict


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


# def _add_overlapping_peaks_single_df(df_regions: pd.DataFrame, peaks_df: pd.DataFrame, r: int) -> pd.DataFrame:
#     """
#     Add column with vector of overlapping peaks.
#     Args:
#         df_regions: pandas DataFrame with regions.
#         peaks_df: pandas DataFrame with peaks.
#         r: radius of the region around the loop center.
#     Returns:
#         pandas DataFrame with added column with vector of overlapping peaks.
#     """
#     main_df = copy.deepcopy(df_regions)
#     main_df['open_regions'] = [[]]*len(main_df)
#     # group by chromosome
#     main_df_grouped = main_df.groupby('chr')
#     peaks_df_grouped = peaks_df.groupby('chr')

#     for main_chr, main_df_chr in tqdm(main_df_grouped):
#         for peaks_chr, peaks_df_chr in peaks_df_grouped:
#             if main_chr == peaks_chr:
#                     for idx, row in main_df_chr.iterrows():
#                             # TO CHANGE !!! - to overlapping peaks
#                             # overlapping_peaks = peaks_df_chr[(peaks_df_chr['start'] >= row[region] - r)&
#                             #                                     (peaks_df_chr['start'] <= row[region] + r)&
#                             #                                     (peaks_df_chr['end'] >= row[region] - r)&
#                             #                                     (peaks_df_chr['end'] <= row[region] + r)]
#                             overlapping_peaks = peaks_df_chr[~(peaks_df_chr['start'] > row['region'] + r)&
#                                                              ~(peaks_df_chr['end'] < row['region'] - r)]
#                             overlapping_peaks['start'] = overlapping_peaks.apply(lambda x: row['region'] - r if x['start'] < row['region'] - r else x['start'], axis=1)
#                             overlapping_peaks['end'] = overlapping_peaks.apply(lambda x: row['region'] + r if x['end'] > row['region'] + r else x['start'], axis=1)
#                             overlapping_peaks['peak'] = tuple(zip(overlapping_peaks['start'], overlapping_peaks['end']))

#                             if len(list(overlapping_peaks['peak'])) > 0:
#                                 main_df.loc[idx, 'open_regions'] = [list(overlapping_peaks['peak'])]
                            
#     return main_df


# def _add_overlapping_peaks(regions_df: pd.DataFrame,
#                           peaks_dfs_dict: Dict[str, pd.DataFrame],
#                           r: int) -> pd.DataFrame:
#     """
#     Add column with vector of overlapping peaks of open chromatin for each dataframe from the main_dfs_dict dictionary.
#     Args:
#         regions_df: pandas DataFrame with regions.
#         peaks_dfs_dict: dictionary with cell types as keys and pandas DataFrames with open chromatin peaks as values.
#         r: radius of the region around the loop center.
#     Returns:
#         dictionary with cell types as keys and pandas DataFrames with regions of anchors and overlapping peaks of open chromatin as values.
#     """
#     print('Adding overlapping peaks of open chromatin...')
#     df_open_chromatin_dict = {}
#     for cell_type, peaks_df in peaks_dfs_dict.items():
#         print(f'...for {cell_type} cell...')
#         # add overlapping peaks
#         df_regions_and_peaks = _add_overlapping_peaks_single_df(regions_df, peaks_df, r)
#         df_open_chromatin_dict[cell_type] = df_regions_and_peaks
#     print('Done!')

#     return df_open_chromatin_dict


# def _find_motif_single_df(df_open_chromatin: pd.DataFrame, 
#                           motif,
#                           genome: dict):
#     """
#     Add two columns with counts of forward and backward motif occurrences in the each anchor in open chromatin regions.
#     """
#     motif_name = motif.name.upper()
#     motif_pssm = motif.pssm

#     if not f'motif_{motif_name}_f' in df_open_chromatin.columns:
#         df_open_chromatin[f'motif_{motif_name}_f'] = 0
#         df_open_chromatin[f'motif_{motif_name}_b'] = 0

#     chromosoms = df_open_chromatin['chr'].unique()
#     for chr in chromosoms:
#         print(f'...for {chr}...')
#         # read specific chromosome from fasta file
#         sequence = genome[chr]

#         df_one_chr = df_open_chromatin[df_open_chromatin['chr'] == chr]

        
#         for idx, row in df_one_chr.iterrows():
#             regions = row['open_regions']
#             for region in regions:
#                 start, end = region
#                 seq_to_search = sequence[start:end+1]
#                 for position, _ in motif_pssm.search(seq_to_search):
#                     if position >= 0:
#                         df_open_chromatin.loc[idx, f'motif_{motif_name}_f'] += 1
#                 for position, _ in motif_pssm.search(seq_to_search[::-1]):
#                     if position >= 0:
#                         df_open_chromatin.loc[idx, f'motif_{motif_name}_b'] += 1
    
#     return df_open_chromatin



# def find_all_motifs(main_dfs_dict: Dict[str, Callable[[], Any]],
#                     peaks_dfs_dict: Dict[str, pd.DataFrame],
#                     motifs_path: str,
#                     genome_path: str,
#                     r: int) -> pd.DataFrame:
#     """
#     """
#     # if main_dfs_dict values are not DataFrames, load them
#     if not isinstance(list(main_dfs_dict.values())[0], pd.DataFrame):
#         main_dfs_dict = _dict_partitions(main_dfs_dict)

#     df2cols = main_dfs_dict[list(main_dfs_dict.keys())[0]]
#     df_regions = _create_one_region_col(df2cols)
#     df_open_chromatin_dict = _add_overlapping_peaks(df_regions, peaks_dfs_dict, r)

#     genome = _genome2dict(genome_path)

#     print('Finding motifs...')
#     motifs_object = motifs.parse(open(motifs_path), 'jaspar')
#     for motif in motifs_object:
#         for cell_type, df in df_open_chromatin_dict.items():
#             print(f'...for {cell_type} cell...')
#             df_open_chromatin_dict[cell_type] = _find_motif_single_df(df, motif, genome)
#     print('Done!')

#     return df_open_chromatin_dict

    