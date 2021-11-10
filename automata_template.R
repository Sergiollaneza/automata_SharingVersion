
# Depedencies -------------------------------------------------------------

# Check if dependencies are installed
mirror = getOption("repos")
mirror["CRAN"] = "https://www.stats.bris.ac.uk/R"
options(repos = mirror)

# CRAN packages
list_packages_CRAN = c("BiocManager, tidyverse, stringr, plyr")
new_packages_CRAN = list_packages_CRAN[!(list_packages_CRAN %in% installed.packages()[,"Package"])]
if(length(new_packages_CRAN)) install.packages(new_packages_CRAN)


list_packages_Bioconductor <- c("TCGAbiolinks")
new_packages_Bioconductor <- list_packages_Bioconductor[!(list_packages_Bioconductor %in% installed.packages()[,"Package"])]
if(length(new_packages_Bioconductor)) library(BiocManager)
if(length(new_packages_Bioconductor)) BiocManager::install(new_packages_Bioconductor, version = 3.8)

# Install the newest version of automata
library(devtools)
install_github("UEA-Cancer-Genetics-Lab/Automata",
               auth_token = "daba61a8aaabf2ec7a0563e1052ea4d4dc740cda",
               ref = "Workshop")

# Load lautomata
library(automata)
library(here)

# Actual code -------------------------------------------------------------

# Choose parameters

# Write here the name of the project you are interested between quoting marks as in the example
proj <- "TCGA-PRAD"

# Write here TRUE if you want to use microarray data or FALSE if you want RNA-seq counts
micro <- FALSE

# Write here an ID you want to assign to your analysis so they are easily distinguisible
autoId <- "templateOutcome"


###### DO NOT CHANGE ANYTHING BELOW THIS LINE ##################################
###############################################################################
###############################################################################


# 1. BARCODES -----

# Get barcodes for the desired option and extract the barcodes
barcodes_df <- getBarcodes(project = proj, microarray = micro, automataID = autoId)
samples <- barcodes_df$Sample
participants <- barcodes_df$Participant

# # 2. DOWNLOAD -----

# Get expression data for the desired option
downloadExpression(project = proj, barcodes = samples, microarray = micro, automataID = autoId)

# Get clinical data for the desired option
downloadClinical(project = proj, barcodes = participants, automataID = autoId)


# Get CNV data for the desired option
downloadCNV(project = proj, barcodes = samples, automataID = autoId)


# Get SNV data for the desired option
downloadSNV(project = proj, barcodes = samples, automataID = autoId)


# Get methylation data for the desired option
downloadMeth(project = proj, barcodes = samples, automataID = autoId)


# 3. PREPROCESS-----

# Removes duplicated samples
cleanMatrix(proj, autoId, micro)

# Normalise transcriptome expression
normaliseExpression(proj, autoId, micro)

# Picks the top 500 genes for LPD
topGenes(proj, autoId, microarray = micro)

# 3. Runs LPD-----

runLPD_A(proj, autoId)

n_jobs <- length(system(paste0("bjobs | grep ", substr(autoId, nchar(autoId) - 9, nchar(autoId))),
                        intern = TRUE))

while(n_jobs > 0){
  message("LPD-A is not finished. Trying again in 30 min...")
  Sys.sleep(1300)
  n_jobs <- length(system(paste0("bjobs | grep ", substr(autoId, nchar(autoId) - 9, nchar(autoId))),
                          intern = TRUE))
}

# Calculate parameters
bestParameters <- estimateParameters(proj, autoId)

# Run LPD-B
runLPD_B(proj, autoId)

n_jobs <- length(system(paste0("bjobs | grep ", substr(autoId, nchar(autoId) - 9, nchar(autoId))),
                        intern = TRUE))

while(n_jobs > 0){
  message("LPD-B is not finished. Trying again in 30 min...")
  Sys.sleep(1300)
  n_jobs <- length(system(paste0("bjobs | grep ", substr(autoId, nchar(autoId) - 9, nchar(autoId))),
                          intern = TRUE))
}

# 4. Postprocess ----
postprocess(proj, autoId)

# 5. Differential analysis ---

# Batch effect
batchEffect(proj, autoId)

# Clinical
clinicalAnalysis(proj, autoId)

# Expression
diffExpGenes(proj, autoId)

# Methylation
methylationAnalysis(proj, autoId)

# CNV
cnvAnalysis(proj, autoId)

# SNV
snvAnalysis(proj, autoId)

# Somatic signatures
somaticSignaturesAnalysis(proj, autoId)

# Dendrogram comparison with LPD
runDendrogram(proj, autoId)

# Genes overlap
createIntersection(proj, autoId)

# 6. REPORT ----
generateReport(proj, autoId, micro)


