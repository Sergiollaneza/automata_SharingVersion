#' Download TCGA data
#'
#' Downloads the expression data from the TCGA and the clinical, methylation, and copy number variation files which are also present in the expression files
#' downloaded. The expression files could be obtained from RNA-seq or microarray experimental strategy but not both at the same time.
#'
#' @param project Name of the TCGA project data to download (eg. TCGA-PRAD or TCGA-BRCA).
#' @param rna.seq  Choose between downloading RNA-seq or Microarray expression data. *TRUE* means RNA-Seq and *FALSE* means Microarray. Defauls to *TRUE*.
#' @param clinical Choose to download or not the clinical data available for the expression files. Defaults to *TRUE*.
#' @param methylation Choose to download or not the metylation data available for the expression files. Defaults to *TRUE*.
#' @param cnv Choose to download or not the copy number variation data available for the expression files. Defaults to *TRUE*.
#' @param save.files Choose to save or not the individual files. Defaults to *TRUE*. If set as FALSE, it will create Robjects with all the data.
#' @param control Choose to download tumour and control samples or only tumour samples. Defatuls to *TRUE*.
#' @export
#' @examples
#' to_do()


# Summary -----------------------------------------------------------------


# This function is part of the Automata package.The whole function is based on the package TCGAbiolinks, therefore an understading of the package would be
# useful for a fully comprehension of how the function works.
# The function is divided into 4 sections, one per each type of file, but the workflow is almost the same for every step. A query with the desired filters
# is used to check which files are available for those requisites. Then these files are downloaded and its content is processed into a data frame to work with it
# in downstream analyses. The main difference of this process between expression and the others types of data, is that for the others, all the found files that are
# not present in the expression list of files, are deprecated. The reason for this is that the main objective of the package is to study gene expression.

# Input: Project name
# Output: Two folders: GDCdata contains the individual files if they are not deleted. TCGA_Results contains CSV files with all the data from the files.
#
# Author: Sergio Llaneza-Lago
# Date: 5th June 2019

download.TCGA.data <- function(project, rna.seq = TRUE, expression = TRUE, clinical = TRUE, methylation = TRUE, cnv = TRUE, snv = TRUE, save.files = TRUE, control = TRUE){

  # Compatiblity check ------------------------------------------------------------

  # Check the version of R
  version = getRversion()
  version = unlist(version)
  if(version[1] < 3){
    stop("Your version of R seems to be too old. Please, use the version 3.5.0 or later", call. = FALSE)
  } else if (version[1] == 3 & version[2] < 5){
    stop("Your version of R seems to be too old. Please, use the version 3.5.0 or later", call. = FALSE)
  }

  # Check that the project argument is correct

  if (rna.seq == TRUE){
    list_of_projects <- c("TCGA-BRCA", 'TCGA-GBM', 'TCGA-OV', 'TCGA-LUAD', 'TCGA-UCEC', 'TCGA-KIRC', 'TCGA-HNSC', 'TCGA-LGG', 'TCGA-THCA', 'TCGA-LUSC',
                         'TCGA-PRAD', 'TCGA-SKCM', 'TCGA-COAD', 'TCGA-STAD', 'TCGA-BLCA', 'TCGA-LIHC', 'TCGA-CESC', 'TCGA-KIRP', 'TCGA-SARC', 'TCGA-LAML',
                         'TCGA-ESCA', 'TCGA-PAAD', 'TCGA-PCPG', 'TCGA-READ', 'TCGA-TGCT', 'TCGA-THYM', 'TCGA-KICH', 'TCGA-ACC', 'TCGA-MESO', 'TCGA-UVM',
                         'TCGA-DLBC', 'TCGA-UCS', 'TCGA-CHOL')

    if((project %in% list_of_projects) == FALSE){
      stop(cat("The project name is not valid, please make sure to use one of the following projects:\n TCGA-BRCA, TCGA-GBM, TCGA-OV, TCGA-LUAD, TCGA-UCEC,
               TCGA-KIRC, TCGA-HNSC, TCGA-LGG, TCGA-THCA, TCGA-LUSC,TCGA-PRAD, TCGA-SKCM, TCGA-COAD, TCGA-STAD, TCGA-BLCA, TCGA-LIHC, TCGA-CESC, TCGA-KIRP, TCGA-SARC,
               TCGA-LAML, TCGA-ESCA, TCGA-PAAD, TCGA-PCPG, TCGA-READ, TCGA-TGCT, TCGA-THYM, TCGA-KICH, TCGA-ACC, TCGA-MESO, TCGA-UVM, TCGA-DLBC, TCGA-UCS, TCGA-CHOL.\n"),
           call.=FALSE)
    }

  } else {
    list_of_projects = c('TCGA-BRCA', 'TCGA-OV', 'TCGA-COAD', 'TCGA-KIRC', 'TCGA-GBM', 'TCGA-LAML', 'TCGA-LUSC', 'TCGA-READ', 'TCGA-UCEC',
                         'TCGA-LUAD', 'TCGA-LGG', 'TCGA-KIRP')

    if((project %in% list_of_projects) == FALSE){
      stop(cat("The project name is not valid or available for micorarray data. Please, switch to RNA-seq data or try one of these projects: 'TCGA-BRCA', 'TCGA-OV',
               'TCGA-COAD', 'TCGA-KIRC', 'TCGA-GBM', 'TCGA-LAML', 'TCGA-LUSC', 'TCGA-READ', 'TCGA-UCEC', 'TCGA-LUAD', 'TCGA-LGG', 'TCGA-KIRP'\n"),
           call.=FALSE)
    }
  }


  # Dependencies ------------------------------------------------------------

  # Configurate a mirror for the packages
  mirror = getOption("repos")
  mirror["CRAN"] = "https://www.stats.bris.ac.uk/R"
  options(repos = mirror)

  # CRAN packages
  list_packages_CRAN = c("tidyr", "here", "BiocManager", "dplyr", "vcfR", "tidyverse", "readr", "purrr",
                         "magrittr", "stringr")
  new_packages_CRAN = list_packages_CRAN[!(list_packages_CRAN %in% installed.packages()[,"Package"])]
  if(length(new_packages_CRAN)) install.packages(new_packages_CRAN)
  library(dplyr)
  library(tidyr)
  library(BiocManager)
  library(here)
  require(vcfR)
  library(tidyverse)
  library(readr)
  library(purrr)
  library(magrittr)
  library(stringr)

  # Biocondcutor packages
  list_packages_Bioconductor = c("TCGAbiolinks", "maftools", "SummarizedExperiment")
  new_packages_Bioconductor = list_packages_Bioconductor[!(list_packages_Bioconductor %in% installed.packages()[,"Package"])]
  if(length(new_packages_Bioconductor)) BiocManager::install(new_packages_Bioconductor, version = 3.8)
  library(TCGAbiolinks)
  library(maftools)
  library(SummarizedExperiment)

  # Creates folder to put all the data
  ifelse(!dir.exists(here("TCGA_results")), dir.create(here("TCGA_results")), FALSE)
  ifelse(!dir.exists(here("TCGA_results", sprintf("%s", project))), dir.create(here("TCGA_results", project)), FALSE)
  ifelse(!dir.exists(here("TCGA_results", sprintf("%s", project), "files")), dir.create(here("TCGA_results", project, "files")), FALSE)
  ifelse(!dir.exists(here("TCGA_results", sprintf("%s", project), "data")), dir.create(here("TCGA_results", project, "data")), FALSE)

  # 1. Download expression data ------------------------------------------------
  cat("\n#\n##\n###\n GETTING EXPRESSION DATA \n###\n##\n#\n")
  if(rna.seq == TRUE){
    expWF = "RNA-seq"

    # 1. 1 RNA-seq Expression data =================================
    # Creates a query that applies different filters to the files and gives back the found files in the TCGA database
    if(control == FALSE){
      expression_query = GDCquery(project = project,
                                  data.category = "Transcriptome Profiling",
                                  data.type = "Gene Expression Quantification",
                                  experimental.strategy = "RNA-Seq",
                                  workflow.type = "HTSeq - Counts",
                                  sample.type = "Primary solid Tumor",
                                  barcode = c("TCGA-*")) # I want all the files that matches that filters
    } else {
      expression_query = GDCquery(project = project,
                                  data.category = "Transcriptome Profiling",
                                  data.type = "Gene Expression Quantification",
                                  experimental.strategy = "RNA-Seq",
                                  workflow.type = "HTSeq - Counts",
                                  barcode = c("TCGA-*")) # I want all the files that matches that filters
    }

    # Generates a table data of the found files
    expression_files = getResults(expression_query)

    # Downloads only if required
    if(expression == TRUE){
      # Downloads the files (it will create a series of subfolders in the working directory automatically)
      cat("RNA-seq expression files are going to be downloaded. To change to microarray please add the argument 'rna.seq=FALSE' to the function.")
      GDCdownload(expression_query,
                  method = "api",
                  files.per.chunk = 70) # Divides the files to download in groups of 70 (if there are enough).

      # Creates a summarised experiment (SE) object which contains the expression matrix (not normalised) among other info (this extra info is not really accurate)
      # and saves the object as a file
      cat("Generating expression matrix.../n")
      SE_expression = GDCprepare(query = expression_query, remove.files.prepared = !(save.files), save = TRUE,
                                save.filename = here("TCGA_results", sprintf("%s", project), "SE_expression_RNA-seq.RData"))

      # Extracts the expression matrix
      expression_matrix = assay(SE_expression)

      # Saves the info in CSV format files
      write.csv(expression_files, here("TCGA_results", project, "files", sprintf("Expression_files_RNA-seq_%s.csv", project)))
      write.csv(expression_matrix, here("TCGA_results", project, "data", sprintf("Expression_matrix_RNA-seq_%s.csv", project)))
    } else{
        cat("Expression files were not downloaded. To change this please add the argument 'expression=TRUE' to the function.")
      }
  }
  if(rna.seq == FALSE){
    expWF = "Microarray"

    # 1. 2 Microarray Expression data =================================

    # Creates a query that applies different filters to the files and gives back the found files in the TCGA database
    if(control == FALSE){
      expression_query = GDCquery(project = project,
                                  data.category = "Gene expression",
                                  data.type = "Gene expression quantification",
                                  experimental.strategy = "Gene expression array",
                                  sample.type = "Primary solid Tumor",
                                  barcode = c("TCGA-*"),
                                  legacy = TRUE)
    }  else{
      expression_query = GDCquery(project = project,
                                  data.category = "Gene expression",
                                  data.type = "Gene expression quantification",
                                  experimental.strategy = "Gene expression array",
                                  barcode = c("TCGA-*"),
                                  legacy = TRUE)
    }

    # Generates a table data of the found files
    expression_files = getResults(expression_query)

    if(expression == TRUE){
      # Downloads the files (it will create a series of subfolders in the working directory automatically)
      cat("Microarray expression files are going to be downloaded. To change to RNA-seq please add the argument 'rna.seq=TRUE' to the function.")
      GDCdownload(expression_query,
                  method = "api",
                  files.per.chunk = 70) # Divides the files to download in groups of 70 (if there are enough).

      # Creates a summarised experiment (SE) object which contains the expression matrix (not normalised) among other info (this extra info is not really accurate)
      # and saves the object as a file
      SE_expression = GDCprepare(query = expression_query, remove.files.prepared = !(save.files), save = TRUE,
                                 save.filename = here("TCGA_results", sprintf("%s", project), "SE_expression_Microarray.RData"))

      # Extracts the expression matrix
      expression_matrix = assay(SE_expression)

      # Saves the info in CSV format files
      # write.csv(expression_files, here("TCGA_results", project, sprintf("Expression_files_Microarray_%s.csv", project))) # Gives some problems with microarray
      write.csv(expression_matrix, here("TCGA_results", project, "data", sprintf("Expression_matrix_Microarray_%s.csv", project)))
    } else {
      cat("Expression files were not downloaded. To change this please add the argument 'expression=TRUE' to the function.")
    }
  }


  # 2. Download clinical data ------------------------------------------------
  if(clinical == TRUE){
    cat("\n#\n##\n###\n DOWNLOAD CLINICAL DATA \n###\n##\n#\n")

    # Creates a query that applies different filters to the files and gives back the found files in the TCGA database
    clinical_query = GDCquery(project = project,
                              data.category = "Clinical",
                              file.type = "xml",
                              barcode = c("TCGA-*"))

    # Generate table with the results and makes sure to have one row for each case (not multiple cases in one row) -> Cases barcodes are in column 8
    clinical_files = getResults(clinical_query)
    clinical_files = separate_rows(clinical_files, 8, sep = ",")

    # Remove the rows which barcode is not present in the barcode of the expression files
    # The barcodes have the 12 first characters exclusive of the patient, and from there, its different per sample (expression, methylation...)
    clinical_matched = clinical_files[which(substr(clinical_files[,8], 0, 12) %in% substr(expression_files[,8], 0, 12)),]

    # Repeats the query but using the barcodes which are overlapped
    clinical_query = GDCquery(project = project,
                              data.category = "Clinical",
                              file.type = "xml",
                              barcode = clinical_matched$cases)

    # Generates a table data of the found files
    clinical_files = getResults(clinical_query)
    clinical_files = separate_rows(clinical_files, 8, sep = ",")

    # Downloads the clinical data
    GDCdownload(clinical_query,
                method = "api",
                files.per.chunk = 70)

    # Creates a data frame with all the patient data and another one with all the sample data (does not give the option to remove the files)
    prep_clinical = GDCprepare_clinic(clinical_query, clinical.info = "patient")
    prep_clinical = unique(prep_clinical)

    # Saves everything as CSV format files
    write.csv(clinical_files, here("TCGA_results", project, "files", sprintf("Clinical_files_%s_%s.csv", project, expWF)))
    write.csv(prep_clinical, here("TCGA_results", project, "data", sprintf("Clinical_data_%s_%s.csv", project, expWF)))

  }

  # 3. Download methylation data ------------------------------------------------
  if(methylation == TRUE){
    cat("\n#\n##\n###\n DOWNLOAD METHYLATION DATA \n###\n##\n#\n")

    # Creates a query that applies different filters to the files and gives back the found files in the TCGA database
    meth_query = GDCquery(project = project,
                          data.category = "DNA Methylation",
                          barcode = "TCGA-*")

    # Generate table with the results and makes sure to have one row for each case (not multiple cases in one row) -> Cases barcodes are in column 8
    meth_files = getResults(meth_query)
    meth_files = separate_rows(meth_files, 8, sep = ",")

    # Remove the rows which barcode is not present in the barcode of the expression files
    # The barcodes have the 12 first characters exclusive of the patient. From 12 to 19, the characters are in common if they are from the same exact sample.
    meth_matched = meth_files[which(unique(substr(meth_files[,8], 0, 19)) %in% substr(expression_files[,8], 0, 19)),]

    # Repeats the query but using the barcodes which are overlapped
    meth_query = GDCquery(project = project,
                          data.category = "DNA Methylation",
                          barcode = meth_matched$cases)

    # Generates a table data of the found files
    meth_files = getResults(meth_query)
    meth_files = separate_rows(meth_files, 8, sep = ",")

    # Downloads the methylation data
    GDCdownload(meth_query,
                method = "api",
                files.per.chunk = 25) # Methylation files are usually big, so the chunks have to be smaller or they will crash

    # Creates a data frame with all the methylation data from all the patients.
    prep_meth = prepare.meth(meth_query, meth_files)

    # Saves everything as CSV format files
    write.csv(meth_files, here("TCGA_results", project, "files", sprintf("Methylation_files_%s_%s.csv", project, expWF)))
    write.csv(prep_meth, here("TCGA_results", project, "data", sprintf("Methylation_data_%s_%s.csv", project, expWF)))
  }

  # 4. Download copy number variation (CNV) data ------------------------------------------------
  if(cnv == TRUE){
    cat("\n#\n##\n###\n DOWNLOAD CNV DATA \n###\n##\n#\n")

    # Creates a query that applies different filters to the files and gives back the found files in the TCGA database
    copy_query = GDCquery(project = project,
                          data.category = "Copy Number Variation",
                          data.type = "Copy Number Segment",
                          barcode = c("TCGA-*"))

    # Generate table with the results and makes sure to have one row for each case (not multiple cases in one row) -> Cases barcodes are in column 8
    copy_files = getResults(copy_query)
    copy_files = separate_rows(copy_files, 8, sep = ",")

    # Remove the rows which barcode is not present in the barcode of the expression files
    # The barcodes have the 12 first characters exclusive of the patient. From 12 to 19, the characters are in common if they are from the same exact sample.
    copy_matched = copy_files[which(unique(substr(copy_files[,8], 0, 19)) %in% substr(expression_files[,8], 0, 19)),]

    # Repeats the query but using the barcodes which are overlapped
    copy_query = GDCquery(project = project,
                          data.category = "Copy Number Variation",
                          data.type = "Copy Number Segment",
                          barcode = copy_matched$cases)

    # Generates a table data of the found files
    copy_files = getResults(copy_query)
    copy_files = separate_rows(copy_files, 8, sep = ",")

    # Downloads the copy number variation data
    GDCdownload(copy_query,
                method = "api",
                files.per.chunk = 25)

    # Creates a data frame with all the CNV data from all the patients.
    prep_copy = GDCprepare(query = copy_query, remove.files.prepared = !(save.files), save = TRUE,
                           save.filename = here("TCGA_results", sprintf("%s", project), "SE_CNV.RData"))

    # Saves everything as CSV format files
    write.csv(copy_files, here("TCGA_results", project, "files", sprintf("cnv_files_%s_%s.csv", project, expWF)))
    write.csv(prep_copy, here("TCGA_results", project, "data", sprintf("cnv_data_%s_%s.csv", project, expWF)))
  }

  # 5. Download SNV data (SNPs + indels) ------------------------------------------------
  if(snv == TRUE){
    cat("\n#\n##\n###\n DOWNLOAD SNV DATA \n###\n##\n#\n")
  }

  # SNV data comes in four different formats. The consensus of the four of them will be selected.
  snv_types <- c("MuSE", "MuTect2", "SomaticSniper", "VarScan2")
  for(i in 1:length(snv_types)){

  # Creates a query that applies different filters to the files and gives back the found files in the TCGA database
  snv_query = GDCquery(project = project,
                       data.category = "Simple Nucleotide Variation",
                       experimental.strategy = "WXS",
                       workflow.type = snv_types[i],
                       access = "controlled",
                       barcode = c("TCGA-*"))

  # Generate table with the results and makes sure to have one row for each case (not multiple cases in one row) -> Cases barcodes are in column 8
  snv_files = getResults(snv_query)
  snv_files <- separate_rows(snv_files, 8, sep = ",")

  # Remove the rows which barcode is not present in the barcode of the expression files
  # The barcodes have the 12 first characters exclusive of the patient. From 12 to 19, the characters are in common if they are from the same exact sample.
  snv_matched = snv_files[which(unique(substr(snv_files[,8], 0, 19)) %in% substr(expression_files[,8], 0, 19)),]

  # Repeats the query but using the barcodes which are overlapped
  snv_query = GDCquery(project = project,
                       data.category = "Simple Nucleotide Variation",
                       experimental.strategy = "WXS",
                       workflow.type = snv_types[i],
                       access = "controlled",
                       barcode = snv_matched$cases)

  # Generates a table data of the found files
  snv_results = getResults(snv_query)
  snv_files = separate_rows(snv_results, 8, sep = ",")

  # Donwloads the SNV data
  GDCdownload(snv_query,
              token.file = token_file,
              method = "client",
              files.per.chunk = 25)

  prep_snv = prepare.snv(snv_query, snv_files)

  # Saves everything as CSV format files
  write.csv(snv_files, here("TCGA_results", project, "files", sprintf("snv_%s_files_%s.csv", snv_types[i], project)))
  write.csv(prep_snv, here("TCGA_results", project, "data", sprintf("snv_%s_data_%s.csv", snv_types[i], project)))
  }

  consensus <- create.consensus(project, snv_types)
  write_csv(consensus, here("TCGA_results", project, "data", sprintf("snv_consensus_%s_%s.csv", project, expWF)))




  # 6. A bit of cleaning ----------------------------------------------------

  if(file.exists(here("MANIFEST.txt")) == T){
    file.remove("MANIFEST.txt")
  }
}
