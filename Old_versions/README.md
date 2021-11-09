# Automata 1.1 - *Now with extra of SNVs*

Automata (**Automa**tic **T**CG**A** data download and processing package) is a R package for downloading and processing TCGA data in the UEA's HPC. So far, the package is able to:

1. Download data. The package is able to download expression files, clinical files, methylation files, copy number variation (CNV) files, and single nucleotide variant (SNV) data. It must be noticed that the cases present in the expression pool of data are used as an index for the package to know which for which cases has to download the other data.

2. Merge the data from different samples into one data frame stored as a CSV file for each type of data.

**These functions are temporarily disabled**:

3. Normalise the expression data by *Variance Stabilizing Transformation* and *log2fold*. Only the 500 genes with the highest variance are selected for the following steps.

4. Launch the first stage of LPD to calculate which combination of number of processes and sigma value works better.

## Getting Started
These instructions will get Automata installed in your system.

### Prerequisites

1. R installed in your system. You can download R from [here](https://www.r-project.org/).

2. A GitHub account with access to the UEA Cancer Genomics Lab's repositories. You can access the team GitHub [here](https://github.com/UEA-Cancer-Genetics-Lab). *Let's ignore the fact that you need to have access to read this anyway.*

3. A GitHub token associated with the account with acess to the repository. You can find a guide to help you to create it here: *just joking I didn't have time to write it yet*

4. The CRAN package **devtools** installed in your system. To install and load **devtools** package, run this in your R:
```
install.packages("devtools")
library(devtools)
```

5. If SNV data is needed  an account in the GDC is required. You can find more info [here](https://gdc.cancer.gov/access-data/obtaining-access-controlled-data). If you already have access, down load the token from the TCGA database.


### Installing

If **devtools** is not loaded, load it.
```
library(devtools)
```
Run the following command to install the package, you will need to add your own personal GitHub token.
```
install_github("UEA-Cancer-Genetics-Lab/Automata", auth_token = "replace_this_with_your_token")
```

Well done! Automata is now succesfully installed in your system.

### Running Automata

Right now the available functions are:

1. Download and merge TCGA data. You can find more info here *in progress*.
```
download.TCGA.data(project, rna.seq = TRUE, expression = TRUE, clinical = TRUE, methylation = TRUE, cnv = TRUE, snv = TRUE, save.files = TRUE)
```
2. That's it.


## Authors

* **Sergio Llaneza-Lago** - *Developer* - [GitHub profile](https://github.com/Sergiollaneza).


## Acknowledgments

* The [developer team](http://bioconductor.org/packages/release/bioc/html/TCGAbiolinks.html) of **TCGAbiolinks**.
* [ChrisEllis93](https://github.com/ChrisEllis93) for testing my package.
