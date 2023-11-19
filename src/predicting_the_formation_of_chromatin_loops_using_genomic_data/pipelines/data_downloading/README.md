# Pipeline data_downloading

## Overview

Pipeline used to download data. The url of the data to be downloaded should be provided in the configuration file: conf/base/parameters/data_downloading.yml

## Configuration

The configuration file is an yml file and is divided into 2 sections: "to_download" and "dm3todm6_mapping". The first section is responsible for downloading data from the given urls. The second section is responsible for mapping the specified  files for the fruit fly from the reference genome version dm3 to dm6. Each section has its own structure:

1. "to_download" - subsequent levels in the structure of this section correspond to the directories being created
    - organism (created directory: data/01_raw/organism)
        - type of data  (created directory: data/01_raw/organism/type_of_data)
            - specific file (name of downloaded file: data/01_raw/organism/type_of_data/specific_file)
                - url to download the specific file

2. "dm3todm6_mapping" - subsequent levels indicate the type of data and the specific file from that data type to be mapped to the dm6 reference genome
    - type of data
        - specific file


## Example

Running only the data download pipeline: 

    kedro run -p data_downloading 