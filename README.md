# Automata Sharing Version

This a version of my package Automata with that does not contain confidential data, therefore, some of the functions and data have been removed. This version is just for display and has been not properly tested for professionally working with it. Even so, please do not share any of the code since this is still a work in progress.

Automata is a package for the downloading and processing of data from The Cancer Genome Atlas database and right now its use is restricted for the Cancer Genetics deparment at the University of East Anglia.

The workflow of the package is divided into six steps: (1) Data downloading,  (2) Data pre-processing, (3) Running LPD, (4) LPD Postprocessing, (5) Differential analysis and (6) Report. We strongly recommend to follow the workflow from beggining to end (you can find a script for that [right here](https://github.com/UEA-Cancer-Genetics-Lab)), but the package allows using just a selection of functions or using custom data. Please, notice that the LPD-related functions are meant to be run on the HPC cluster environment of the UEA.


## Data downloading

Automata downloads the data from The Cancer Genome Atlas an can be used for transcriptome (RNA-seq or microarray), clinical, methylation, SNV, and CNV data.

The TCGA is structured into projects and each contains the data belonging to a major cancer type e.g. prostate adenocarcinoma data is on TGGA-PRAD. Projects are made up of cases which represent a unique patient and each case contains multiple files such as transcriptome or clinical features among others. 

Therefore, to download the data, it is required to have the barcodes representing each patient. For this we can use the `getBarcodes` function and then use its output on the corresponding downloading function. Here is an example of how we can download all the data for prostate cancer.









