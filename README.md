# Automata Sharing Version

This a version of my package Automata with that does not contain confidential data, therefore, some of the functions and data have been removed. This version is just for display and has been not properly tested for professionally working with it. Even so, please do not share any of the code since this is still a work in progress.

Automata is a package for the downloading and processing of data from The Cancer Genome Atlas database and right now its use is restricted for the Cancer Genetics deparment at the University of East Anglia.

The workflow of the package is divided into six steps: (1) Data downloading,  (2) Data pre-processing, (3) Running LPD, (4) LPD Postprocessing, (5) Differential analysis and (6) Report. We strongly recommend to follow the workflow from beggining to end (you can find a script for that [right here](https://github.com/UEA-Cancer-Genetics-Lab)), but the package allows using just a selection of functions or using custom data. Please, notice that the LPD-related functions are meant to be run on the HPC cluster environment of the UEA.


## 1. Data downloading

Automata downloads the data from The Cancer Genome Atlas an can be used for transcriptome (RNA-seq or microarray), clinical, methylation, SNV, and CNV data.

The TCGA is structured into projects and each contains the data belonging to a major cancer type e.g. prostate adenocarcinoma data is on TGGA-PRAD. Projects are made up of cases which represent a unique patient and each case contains multiple files such as transcriptome or clinical features among others. 

Therefore, to download the data, it is required to have the barcodes representing each patient. For this we can use the `getBarcodes` function and then use its output on the corresponding downloading function. Here is an example of how we can download all the data for prostate cancer selecting RNA-seq as preferred transcriptome format.

``` {r, eval = FALSE}

library(automata)

# Fetches the barcodes, sample codes and cases code from the TCGA
barcodes_df <- getBarcodes(project = "TCGA-PRAD", microarray = FALSE, automataID = "example")

# Extracts samples and participants as vectors
samples <- barcodes_df$Sample
participants <- barcodes_df$Participant

# Get expression data for prostate cancer
downloadExpression(project = "TCGA-PRAD", barcodes = samples, microarray = micro, automataID = "example")

# Get clinical data for prostate cancer
downloadClinical(project = "TCGA-PRAD", barcodes = participants, automataID = "example")


# Get CNV data for prostate cancer
downloadCNV(project = "TCGA-PRAD", barcodes = samples, automataID = "example")


# Get SNV data for prostate cancer
downloadSNV(project = "TCGA-PRAD", barcodes = samples, automataID = "example")


# Get methylation data for prostate cancer
downloadMeth(project = "TCGA-PRAD", barcodes = samples, automataID = "example")

```

## 2. Preprocessing of the transcriptome data

TCGA transcriptome data is available as RNA-seq counts or microarray intensity quantification. Both methodologies require to apply normalisation in order to be comparable across samples and they need to be adapted for LPD. 

The function `cleanMatrix` removes duplicate transcriptome samples and we strongly recommend using it if your goal is to apply LPD. Duplicated samples can originate errors downstream when the transcriptome data is matched with other data such as methylation levels or clinical features. Then `normaliseExpression` normalises RNA-seq counts by using Variance Stabilizacion Transformation and applying logarithms, while microarray intensities are normalised by adding the most negative to remove negative values. Finally, `topGenes` pick the 500 genes with the highest variance across all samples as it is defined in the LPD protocol.

Following the previous example, here we apply the preprocessing functions to the expression matrix we obtained previously. 

## 3. Running LPD

*The functions described on this step were developed to be used on the HPC of the University of East Anglia where LPD is installed.*

LPD is a Bayesian unsupervised clustering algorithm that is able to classify expression samples into soft clusters taking into account the high genetic heterogeneity of cancer samples.

LPD is run in two stages: The first one (Stage A) calculates the fitness of several hidden parameters combinations by using log-likelihood and picks the three best candidates. The second stage (Stage B) runs each combination 100 times, but due to job queue limits need to be run for several days.

Once LPD is done, the output is a matrix with a row for each sample and a column for each detected number of subtypes (or clusters) that contains a value from 0 to 1 representing how much each sample is represented by each subtype.

## 4. Postprocessing

After running LPD, the script produces 300 runs (100 for each optimal candidate). The function `postprocess` compares all the runs and picks the best one according to the internal correlations for each run. Then, it assigns each sample to a cluster, and it is assumed that each cluster represents a unique subtype or molecular signature that explains the genetic variability across all samples.

To execute the postprocess function as we did before, just do:

## 5. Differential analysis

The differential analysis consists on comparing each subtype to the others with the purpose of identifying the unique traits of all the subtypes. Here we can use all the plethora of data available in the TCGA that was downloaded in the first step.

In this example, we run all possible differential analysis for prostate cancer.

## 6. Report

The `generateReport` function uses Rmarkdown and Flexdashboard to create an static HTML file with a dashboard that displays a summary of the analysis outcome. 









