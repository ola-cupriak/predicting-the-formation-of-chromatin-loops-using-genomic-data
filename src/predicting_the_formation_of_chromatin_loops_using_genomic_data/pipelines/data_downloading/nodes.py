import subprocess
import os
import yaml
    
    
def _create_dir(dir_path):
    """
    Creates a directory if it does not exist.
    ::param dir_path: path to the directory
    """
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)


def download_data(datasources_dict: dict):
    """
    Downloads the data from the urls in the datasources_dict 
    to created experiments directories.
    Creates a yml file in which the names of downloaded files correspond to cell types.
    ::param datasources_dict: dictionary with the urls to download the data
    """
    names_dict = {}
    for experiment, cell_dict in datasources_dict.items():
        names_dict[experiment] = {}
        # create directory for the experiment
        dir_experiment = 'data/01_raw/'+experiment
        _create_dir(dir_experiment)
        # change directory to the experiment directory
        os.chdir(dir_experiment)
        files_before = set(os.listdir())
        gunzip = False
        for cell_type, url in cell_dict.items():
            # download the data
            subprocess.run(['wget', url])
            files_after = set(os.listdir())
            file_name = list(files_after - files_before)[0]
            if file_name.endswith(".gz"):
                gunzip = True
                file_name = file_name[:-3]
            names_dict[experiment][cell_type] = file_name
            files_before = files_after
        # gunzip
        if gunzip:
            subprocess.run(f'gunzip {" ".join(list(files_after))}', shell=True)
        # change directory to the experiment directory
        os.chdir('../../..')
    yaml_string=yaml.dump(names_dict)
    with open('data/01_raw/cells2names.yml', 'w') as f:
        f.write(yaml_string)