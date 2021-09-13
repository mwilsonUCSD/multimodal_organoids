# Introduction
This code was used to analyze electrophysiology and imaging data from cortical organoids implanted in NOD/SCID mice brains and surrounding cortex.

# Publication
TBD

# Usage Guide

## Requirements
This code was written and tested using MATLAB v2019b. Two external toolboxes are required to run the code:
- Chronux toolbox (http://chronux.org/)
- Circular Statistics toolbox (https://www.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics)

To use the external toolboxes, download the toolbox and save it to the main folder then add the necessary functions to your MATLAB file path.

## Licence
TBD

## Instructions
As a first step, convert the raw electrophysiology data (recorded with an Intan RHD2000 recording amplifier) by running the function 'read_Intan_RHD2000_file.m' with the file name as the function input. After this step is complete, the data may be saved as a .mat file for quicker future reference.

There are three main MATLAB files included in this folder: one for each type of analysis (LFP, MUA, and Anesthesia) and their corresponding figures. The sections within each file are meant to be run sequentially and call on helper functions included in the folder 'helper functions.' For each file, the first section asks you to enter a data file name (see previous paragraph on converting raw data to .mat) then uploads the corresponding file (a single recording from a single day) to be analyzed in subsequent sections.
