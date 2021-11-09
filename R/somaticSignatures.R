
#' Analysis of somatic signatures
#'
#' Does an analysis of somatic signatures
#' @param project
#' @param automataID
#' @import BSgenome.Hsapiens.UCSC.hg38
#' @import magrittr
#' @import ggplot2
#' @importFrom readr read_csv
#' @importFrom dplyr select mutate filter bind_rows summarise_if group_by ungroup
#' @importFrom deconstructSigs mut.to.sigs.input
#' @importFrom MutationalPatterns cos_sim_matrix plot_cosine_heatmap plot_contribution plot_contribution_heatmap fit_to_signatures plot_96_profile
#' @export

somaticSignaturesAnalysis <- function(project, automataID){

  message("\n\n")
  message("-------------------------------------------------------")
  message("------------------SOMATIC SIGNATURES-------------------")
  message("-------------------------------------------------------")
  message("\n\n")

  # Read the assigned samples
  gammaValues <- read_csv(here::here("TCGA_results", project, automataID, "lpd_parameters",
                                     paste0("Sample_assigned_", project, "_", automataID, ".csv")))

  processes <- unique(gammaValues$Max_LPD)

  # Read SNV data
  mutations <- read_csv(here::here("TCGA_results", project, automataID, "TCGA_data",
                                   paste0("SNV_data_", project, "_", automataID, ".csv"))) %>%
    dplyr::select(Tumor_Sample_Barcode, Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2) %>%
    mutate(Tumor_Sample_Barcode = substr(Tumor_Sample_Barcode, 1, 16))

  # Loads Homo Sapiens genome (hg38)
  homoGenome <- BSgenome.Hsapiens.UCSC.hg38

  # Creates matrix dataset for all the samples ------------------------
  mutationsTibble <- tibble()
  tempRowName <- list()
  i = 1

  for(process in processes){
    message("o Starting process ", process)

    # Select samples
    selectedSamples <- dplyr::filter(gammaValues, Max_LPD == process)$Sample

    # Skips groups of only control samples
    if(all(substr(selectedSamples, nchar(selectedSamples[1]) - 2, nchar(selectedSamples[1])) == "11A")){
      next
    }

    mutationsSelected <- dplyr::filter(mutations, Tumor_Sample_Barcode %in% selectedSamples) %>%
      mutate(Tumor_Sample_Barcode = paste0(process, "_", Tumor_Sample_Barcode)) %>%
      as.data.frame()

    if(nrow(mutationsSelected) == 0){
      next
    }

    # Creates the signature matrix
    mutationsMatrix <- mut.to.sigs.input(mut.ref = mutationsSelected,
                                         sample.id = "Tumor_Sample_Barcode",
                                         chr = "Chromosome",
                                         pos = "Start_Position",
                                         ref = "Reference_Allele",
                                         alt = "Tumor_Seq_Allele2",
                                         bsg = homoGenome)

    # Stores the rowname (it gets deleted in the tibble)
    tempRowName[[i]] <- as.data.frame(rownames(mutationsMatrix))

    # Stores the matrix
    mutationsTibble<- bind_rows(mutationsTibble, mutationsMatrix)

    i = i + 1


  }

  # Creates a column with all the rownames and places it on the tibble rownames
  rowsNames <- bind_rows(tempRowName)
  rownames(mutationsTibble) <- rowsNames$`rownames(mutationsMatrix)`
  signaturesProportion <- as.data.frame(t(mutationsTibble))

  # Exports the matrix
  write.csv(signaturesProportion, here::here("TCGA_results", project, automataID, "mutationalSignatures",
                                             paste0("signaturesProportion_", project, "_", automataID, ".csv")))



  # Uses cosmic to check the results ---------------------------

  # Download mutational signatures from COSMIC
  cosmicURL <- paste("https://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt", sep = "")
  cancerSignatures = read.table(cosmicURL, sep = "\t", header = TRUE)

  # Match the order of the mutation types to the standard one of MutationalPatterns
  new_order = match(row.names(signaturesProportion), cancerSignatures$Somatic.Mutation.Type)

  # Reorder cancer signatures dataframe
  cancerSignatures <- cancerSignatures[as.vector(new_order),]

  # Adds trinucleotide changes names as row.names
  row.names(cancerSignatures) = cancerSignatures$Somatic.Mutation.Type

  # Keep only 96 contributions of the signatures in the matrix
  cancerSignatures <- as.matrix(cancerSignatures[,4:33])

  # Process the signatures proportion to get a group representation instead
  signaturesPropGroup <- signaturesProportion %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    mutate(rowname = substr(rowname, 1, 5)) %>%
    group_by(rowname) %>%
    summarise_if(is.numeric, sum, na.rm = TRUE) %>%
    ungroup() %>%
    remove_rownames() %>%
    column_to_rownames() %>%
    t()

  # Calculates similarity between my mutational profile and the COSMIC one
  cosineSimilarity <- cos_sim_matrix(signaturesPropGroup, cancerSignatures)

  pdf(here::here("TCGA_results", project, automataID, "mutationalSignatures",
                 paste0("cosineSimilarity_", project, "_", automataID, ".pdf")),
      width = 15, height = 9, onefile = FALSE)
  print(plot_cosine_heatmap(cosineSimilarity))
  dev.off()

  # Fit my results with COSMIC signatures
  fitSignatures <- fit_to_signatures(signaturesPropGroup, cancerSignatures)

  # Plots the ones with more than 10
  select <- which(rowSums(fitSignatures$contribution) > 10)

  pdf(here::here("TCGA_results", project, automataID, "mutationalSignatures",
                 paste0("contributionBarplot_", project, "_", automataID, ".pdf")),
      width = 15, height = 9, onefile = FALSE)
  print(plot_contribution(fitSignatures$contribution[select,],
                    cancerSignatures[,select], mode = "relative"))
  dev.off()

  # Plots contribution heatmap
  pdf(here::here("TCGA_results", project, automataID, "mutationalSignatures",
                 paste0("contributionHeatmap_", project, "_", automataID, ".pdf")),
      width = 15, height = 9, onefile = FALSE)
  print(plot_contribution_heatmap(fitSignatures$contribution, method = "complete"))
  dev.off()

  # Compares the reconstructed from the fitting with the original one
  cosMatrix <- as.data.frame(diag(cos_sim_matrix(signaturesPropGroup, fitSignatures$reconstructed)))
  colnames(cosMatrix) <- "Similarity"
  cosMatrix$Group <- rownames(cosMatrix)

  pdf(here::here("TCGA_results", project, automataID, "mutationalSignatures",
                 paste0("ReconstructedSimilarity_", project, "_", automataID, ".pdf")),
      width = 9, height = 7, onefile = FALSE)
  print(ggplot(cosMatrix, aes(x = Group, y = Similarity)) +
    geom_bar(stat = "identity", fill = "skyblue4", width = 0.8) +
    coord_flip(ylim=c(0.8,1)) +
    ylab("Cosine similarity\n original VS reconstructed") +
    xlab("") +
    xlim(rev(levels(factor(cosMatrix$Group)))) +
    theme_bw() +
    theme(panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank()) +
    geom_hline(aes(yintercept=.95)) +
    geom_text(aes(label=round(Similarity, 3)), vjust=0.5, hjust = 1.5, color="white", size=3.5))
  dev.off()

  message("\n\n")
  message("=======================================================")
  message("\n\n")

}
