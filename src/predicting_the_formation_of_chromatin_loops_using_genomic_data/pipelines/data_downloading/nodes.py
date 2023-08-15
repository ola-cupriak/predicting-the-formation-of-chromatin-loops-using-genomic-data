from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
import subprocess
import yaml


def _create_dir(dir_path: str) -> None:
    """
    Creates a directory if it does not exist
    and add .gitkeep file to it.
    Args:
        dir_path: path to the directory
    """
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
    with open(dir_path + "/.gitkeep", "w") as f:
        pass


def download_data(datasources_dict: dict) -> str:
    """
    Downloads the data from the urls in the datasources_dict
    to created experiments directories.
    Creates a yml file with dictionary in the template: {dataset_name: {file_name: cell_type}}}
    Args:
        datasources_dict: dictionary with the urls to download the data
    Returns:
        empty (dummy) string
    """
    to_return = ["None", "None", {}, {}, "None"]
    for organism, data_dict in datasources_dict.items():
        names_dict = {}
        for experiment, cell_dict in data_dict.items():
            names_dict[experiment] = {}
            # create directory for the experiment
            dir_experiment = f"data/01_raw/{organism}/{experiment}"
            _create_dir(dir_experiment)
            # change directory to the experiment directory
            os.chdir(dir_experiment)
            files_before = set(os.listdir())
            gunzip = False
            for cell_type, url in cell_dict.items():
                # download the data
                subprocess.run(["wget", url])
                files_after = set(os.listdir())
                file_name = list(files_after - files_before)[0]
                if file_name.endswith(".gz"):
                    gunzip = True
                    file_name = file_name[:-3]
                names_dict[experiment][file_name] = cell_type
                files_before = files_after

                if (
                    organism == "Homo_sapiens"
                    and experiment == "reference_genomes"
                    and cell_type == "hg19"
                ):
                    to_return[0] = f"data/01_raw/{organism}/{experiment}/{file_name}"
                if (
                    organism == "D_melanogaster"
                    and experiment == "reference_genomes"
                    and cell_type == "dm6"
                ):
                    to_return[1] = f"data/01_raw/{organism}/{experiment}/{file_name}"
                if (
                    organism == "D_melanogaster"
                    and experiment == "reference_genomes"
                    and cell_type == "dm3todm6_chain"
                ):
                    to_return[4] = f"data/01_raw/{organism}/{experiment}/{file_name}"
            # gunzip
            if gunzip:
                subprocess.run(f'gunzip {" ".join(list(files_after))}', shell=True)
            # change directory to the experiment directory
            os.chdir("../../../..")
        yaml_string = yaml.dump(names_dict)
        with open(f"data/01_raw/{organism}/cells2names.yml", mode="w") as f:
            if organism == "D_melanogaster":
                f.write(
                    "HiC_loops_annotations:\n\tlong_and_short_range_loops_D_mel.tsv: neuronal\n"
                )
            f.write(yaml_string)
        if organism == "Homo_sapiens":
            to_return[2] = names_dict
        elif organism == "D_melanogaster":
            to_return[3] = names_dict

    return to_return


def simplify_human_genome_file(path: str) -> list:
    """Get human chromosomes dict from fasta file.
    Args:
        path (str): path to fasta file with human genome.
    Returns:
        list of SeqRecord objects.
    """
    if path == "None":
        return []
    genome = SeqIO.parse(open(path), "fasta")
    chromosomes = []
    for sequence in genome:
        desc = sequence.description
        seq_id = sequence.id
        if seq_id.startswith("NC") and "chromosome" in desc:
            chrom = "chr" + desc.split("chromosome ")[1].split(",")[0]
            chromosomes.append(SeqRecord(sequence.seq, id=chrom, description=""))

    return chromosomes


def simplify_flies_genome_file(path: str) -> list:
    """Get flies chromosomes dict from fasta file.
    Args:
        path (str): path to fasta file with dm6 genome.
    Returns:
        list of SeqRecord objects.
    """
    if path == "None":
        return []
    flies_chromosomes = ["chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX", "chrY"]
    genome = SeqIO.parse(open(path), "fasta")
    chromosomes = []
    for sequence in genome:
        desc = sequence.description
        seq_id = sequence.id
        if desc in flies_chromosomes:
            chrom = desc
            chromosomes.append(SeqRecord(sequence.seq, id=chrom, description=""))

    return chromosomes


def map_dm3_to_dm6(cells2names: dict, dm3todm6_mapping: dict, path_chain: str) -> None:
    """
    Maps bigwig files from dm3 to dm6 using CrossMap.py
    Args:
        cells2names: dictionary with the names of the files
        dm3todm6_mapping: dictionary with the files to map
        path_chain: path to the chain file
    """
    try:
        for experiment, cell_list in dm3todm6_mapping.items():
            if experiment not in cells2names.keys():
                print(f"WARNING! No data for {experiment} to map.")
                continue
            filenames = cells2names[experiment]
            filenames = {v: k for k, v in filenames.items()}
            for cell in cell_list:
                filename = filenames[cell]
                path = f"data/01_raw/D_melanogaster/{experiment}/{filename}"
                subprocess.run(
                    f"CrossMap.py bigwig {path_chain} {path} {path}.mapped", shell=True
                )
                os.remove(path)
                os.rename(f"{path}.mapped.bw", path)
                to_rm = path.split("/")[:-1]
                to_rm = "/".join(to_rm)
                subprocess.run(f"rm {to_rm}/*.bgr", shell=True)
                print(f"{experiment} for cell {cell} mapped successfully.\n")
    except:
        pass
