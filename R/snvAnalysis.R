#' SNV analysis
#'
#' Analyse the SNV data downloaded from the TCGA.
#' @param project Name of the project to analyse according to the TCGA (e.g. TCGA-PRAD)
#' @param automataID Custom ID tag
#' @importFrom readr read_csv
#' @importFrom dplyr select rename
#' @import here
#' @importFrom maftools oncoplot read.maf titv
#' @import ggplot2
#' @import magrittr
#' @importFrom ggpubr ggarrange
#' @importFrom reshape2 melt
#' @export

snvAnalysis <- function(project, automataID){

  message("\n\n")
  message("-------------------------------------------------------")
  message("----------------------SNV ANALYSIS---------------------")
  message("-------------------------------------------------------")
  message("\n\n")

  createTCGAfolder(project, automataID)

  # Import ------------------------------------------------------------------

  # I need the consensus, the gammadata with the groups and the clinical data for the plots
  consensus <- read_csv(here::here("TCGA_results", project, automataID, "TCGA_data",
                                   paste0("SNV_data_", project, "_", automataID, ".csv"))) %>%
    mutate(Tumor_Sample_Barcode = substr(Tumor_Sample_Barcode, 1, 16),
           Tumor_Sample_Barcode = gsub("\\.", "-", Tumor_Sample_Barcode))

  clinical <- read_csv(here::here("TCGA_results", project, automataID, "TCGA_data",
                                  paste0("clinical_data_", project, "_", automataID, ".csv"))) %>%
    dplyr::select(-c(1)) %>%
    dplyr::rename(Tumor_Sample_Barcode = bcr_patient_barcode)

  gammaValues <- read_csv(here::here("TCGA_results", project, automataID, "lpd_parameters",
                                     paste0("Sample_assigned_", project, "_", automataID, ".csv")))

  # Extracts all the processes found for the dataset
  processes <- unique(gammaValues$Max_LPD)


  for(process in processes){

    message("o Analysis of process ", process)


    # Creates a vector for the selected samples and the others
    selectedSamples <- filter(gammaValues, Max_LPD == process)$Sample
    contextSamples <- filter(gammaValues, Max_LPD != process)$Sample

    # Check if there is any sample that is not control
    if(length(grep("-01A", selectedSamples)) == 0) next
    if(nrow(dplyr::filter(consensus, Tumor_Sample_Barcode %in% selectedSamples)) == 0) next

    # Creates a MAF object for each one
    mafSelected <- read.maf(maf = filter(consensus, Tumor_Sample_Barcode %in% selectedSamples),
                            clinicalData = filter(clinical, Tumor_Sample_Barcode %in% substr(selectedSamples, 1, 12)),
                            isTCGA = TRUE)

    mafContext <- read.maf(maf = filter(consensus, !(Tumor_Sample_Barcode %in% selectedSamples)),
                           clinicalData = filter(clinical, !(Tumor_Sample_Barcode %in% substr(selectedSamples, 1, 12))),
                           isTCGA = TRUE)


    # Visualize -------------------------------------------------------------

    # Variant classification ------
    variantClass.selected <- mafSelected@variant.classification.summary %>% mutate(Group = process)
    variantClass.context <- mafContext@variant.classification.summary %>% mutate(Group = "Other_LPD")

    variantClass <- bind_rows(variantClass.selected, variantClass.context) %>% replace(is.na(.), 0) %>%
      dplyr::select(Tumor_Sample_Barcode, Group, everything()) %>%
      dplyr::select(-total, everything()) %>%
      filter(total < sd(total) * 3)

    write_csv(variantClass, here::here("TCGA_results", project, automataID, "snv_analysis",
                                       paste0("variantClass_", project, "_", process, "_", automataID, ".csv")))

    variantClass.plotData <- variantClass %>%
      dplyr::select(-total) %>%
      gather(key = "VariantClass", value = "Proportion", 3:ncol(.))

    plot <- ggplot(data = variantClass.plotData,
                   aes(x = reorder(VariantClass, Proportion), y = Proportion, fill = Group)) +
      coord_flip() +
      theme_minimal() +
      ggtitle(paste("Variant classification for ", process)) +
      theme(axis.title.y = element_blank(),
            axis.title.x = element_blank(),
            axis.text.y = element_text(face = "italic"),,
            plot.title = element_text(hjust = 0.5, face = "bold")) +
      geom_boxplot(width=0.1) +
      scale_fill_manual(values=c("#189cc5", '#005f85'), labels = c(process, "Other LPD"))

    ggsave(here::here("TCGA_results", project, automataID, "snv_analysis",
                      paste0("variantClassification_Plot_", project, "_", process, "_", automataID, ".pdf")),
           plot, width = 7, height = 4, units = "in")


    variantClassification_Plot <- plot
    save(variantClassification_Plot, file = paste0("TCGA_results/", project, "/", automataID,
                                                   "/report_plots/variantClassification_Plot_", project, "_", process, "_", automataID, ".Rdata"))

    # Variant type ------
    variantType.selected <- mafSelected@variant.type.summary %>% mutate(Group = process)
    variantType.context <- mafContext@variant.type.summary %>% mutate(Group = "Other_LPD")

    variantType <- bind_rows(variantType.selected, variantType.context) %>% replace(is.na(.), 0) %>%
      dplyr::select(Tumor_Sample_Barcode, Group, everything()) %>%
      dplyr::select(-total, everything())

    write_csv(variantType, here::here("TCGA_results", project, automataID, "snv_analysis",
                                      paste0("variantType_", project, "_", process, "_", automataID, ".csv")))

    variantType.plotData <- variantType %>%
      dplyr::select(-total) %>%
      filter(SNP < sd(SNP) * 3) %>%
      gather(key = "VariantType", value = "Proportion", 3:ncol(.))

    plot <- ggplot(data = variantType.plotData,
                   aes(x = reorder(VariantType, Proportion), y = Proportion, fill = Group)) +
      coord_flip() +
      ylim(0, 100) +
      geom_boxplot(width=0.1) +
      theme_minimal() +
      ggtitle(paste("Variant Type for ", process, project)) +
      theme(axis.title.y = element_blank(),
            axis.title.x = element_blank(),
            axis.text.y = element_text(face = "italic"),,
            plot.title = element_text(hjust = 0.5, face = "bold")) +
      scale_fill_manual(values=c("#189cc5", '#005f85'), labels = c(process, "Other LPD"))

    ggsave(here::here("TCGA_results", project, automataID, "snv_analysis",
                      paste0("variantType_Plot_", project, "_", process, "_", automataID, ".pdf")),
           plot, width = 7, height = 4, units = "in")

    variantType_Plot <- plot
    save(variantType_Plot, file = paste0("TCGA_results/", project, "/", automataID,
                                         "/report_plots/variantType_Plot_", project, "_", process, "_", automataID, ".Rdata"))


    # Variants per sample  ------
    variantSample.selected <- mafSelected@variants.per.sample %>% mutate(Group = process)
    variantSample.context <- mafContext@variants.per.sample %>% mutate(Group = "Other_LPD")

    variantSample <- bind_rows(variantSample.selected, variantSample.context) %>% replace(is.na(.), 0) %>%
      dplyr::select(Tumor_Sample_Barcode, Group, everything())

    write_csv(variantSample, here::here("TCGA_results", project, automataID, "snv_analysis",
                                        paste0("variantSample_", project, "_", process, "_", automataID, ".csv")))

    variantSample.pvalue <- wilcox.test(formula = Variants ~ Group, data = variantSample)$p.value

    plot <- ggplot(variantSample,
                   aes(x = Group, y = Variants,  fill = Group)) +
      geom_violin(trim = FALSE) +
      ylab("Number of variants per sample") +
      coord_flip() +
      theme_minimal() +
      geom_boxplot(width=0.1) +
      scale_fill_manual(values = c("#189cc5", '#005f85')) +
      ggtitle(paste("Variants per sample for", process, project),
              subtitle = paste0("p-value: ", round(variantSample.pvalue, digits = 4))) +
      theme(axis.title.y = element_blank(),
            plot.title = element_text(hjust = 0.5, face = "bold"),
            legend.position = "none")


    ggsave(here::here("TCGA_results", project, automataID, "snv_analysis",
                      paste0("variantSample_Plot_", project, "_", process, "_", automataID, ".pdf")),
           plot, width = 6, height = 4, units = "in")

    variantSample_Plot <- plot

    save(variantSample_Plot, file = paste0("TCGA_results/", project, "/", automataID,
                                           "/report_plots/variantSample_Plot_", project, "_", process, "_", automataID, ".Rdata"))

    # Type of SNP  ------
    snpType.selected <- titv(maf = mafSelected, plot = FALSE)[[1]] %>% mutate(Group = process)
    snpType.context <- titv(maf = mafContext, plot = FALSE)[[1]] %>% mutate(Group = "Other_LPD")

    snpType <- bind_rows(snpType.selected, snpType.context) %>% replace(is.na(.), 0) %>%
      dplyr::select(Tumor_Sample_Barcode, Group, everything())

    write_csv(snpType, here::here("TCGA_results", project, automataID, "snv_analysis",
                                  paste0("snpType_", project, "_", process, "_", automataID, ".csv")))

    snvType.plotData <- snpType %>%
      dplyr::select(-Tumor_Sample_Barcode) %>%
      melt()

    plot <- ggplot(snvType.plotData, aes(x = variable, y = value, fill = Group)) +
      geom_boxplot(width = 0.5) +
      theme_minimal() +
      theme(legend.position = "none",
            axis.title.y = element_text(hjust = 0.5, size = 15),
            axis.title.x = element_blank()) +
      labs(y = "% of mutation") +
      scale_fill_manual(values = c("#189cc5", '#005f85')) +
      ylim(0, 100) +
      ggtitle(paste("SNP type for", process, project))

    ggsave(here::here("TCGA_results", project, automataID, "snv_analysis",
                      paste0("vsnpType_Plot_", project, "_", process, "_", automataID, ".pdf")),
           plot, width = 6, height = 4, units = "in")

    snpType_Plot <- plot

    save(snpType_Plot, file = paste0("TCGA_results/", project, "/", automataID,
                                     "/report_plots/snpType_Plot_", project, "_", process, "_", automataID, ".Rdata"))


    # Oncoplot ------
    # oncoplot_Plot <- oncoplot(mafSelected, showTumorSampleBarcodes = TRUE)


    pdf(here::here("TCGA_results", project, automataID, "snv_analysis",
                   paste0("oncoplot_Plot_", project, "_", process, "_", automataID, ".pdf")),
        width = 8, height = 10)
    oncoplot(mafSelected, showTumorSampleBarcodes = TRUE)
    dev.off()

    pdf(here::here("TCGA_results", project, automataID, "snv_analysis",
                   paste0("oncoplotContext_Plot_", project, "_", process, "_", automataID, ".pdf")),
        width = 8, height = 10)
    oncoPlotOther_Plot <- oncoplot(mafContext)
    dev.off()

    bitmap(here::here("TCGA_results", project, automataID, "snv_analysis",
                   paste0("oncoplot_Plot_", project, "_", process, "_", automataID, ".png")),
        width = 500, height = 750, res = 100)
    oncoplot(mafSelected, showTumorSampleBarcodes = TRUE)
    dev.off()

    bitmap(here::here("TCGA_results", project, automataID, "snv_analysis",
                    paste0("oncoplotContext_Plot_", project, "_", process, "_", automataID, ".png")),
         width = 500, height = 750, res = 100)
    oncoPlotOther_Plot <- oncoplot(mafContext)
    dev.off()


    # Table most mutated genes ------
    mutatedGenesSelected <- mafSelected@gene.summary

    write_csv(mutatedGenesSelected, here::here("TCGA_results", project, automataID, "snv_analysis",
                                               paste0("mutatedGenes_", project, "_", process, "_", automataID, ".csv")))

    # For the context
    mutatedGenesContext <- mafContext@gene.summary

    write_csv(mutatedGenesContext, here::here("TCGA_results", project, automataID, "snv_analysis",
                                               paste0("mutatedGenesContext_", project, "_", process, "_", automataID, ".csv")))

    # Calculates the proportions
    compareMutations <- mutatedGenesSelected %>%
      left_join(mutatedGenesContext, by= "Hugo_Symbol") %>%
      dplyr::select(Hugo_Symbol, MutatedSamples.x, MutatedSamples.y) %>%
      replace(is.na(.) , 0) %>%
      dplyr::rename(MutatedSamples.Selected = MutatedSamples.x, MutatedSamples.Context = MutatedSamples.y) %>%
      mutate(TotalSamples.Selected = as.numeric(mafSelected@summary$summary[3]),
             TotalSamples.Context = as.numeric(mafContext@summary$summary[3]),
             NotMutated.Selected = TotalSamples.Selected - MutatedSamples.Selected,
             NotMutated.Context = TotalSamples.Context - MutatedSamples.Context) %>%
      dplyr::select(Hugo_Symbol, MutatedSamples.Selected, MutatedSamples.Context, NotMutated.Selected, NotMutated.Context, everything())

    compareMutations$chisq <- apply(compareMutations, 1, function(x){
      tbl = matrix(as.numeric(x[2:5]), ncol = 2, byrow = TRUE)
      chisq.test(tbl)$p.value
    })

    compareMutations$padj <- p.adjust(compareMutations$chisq, method = "BH")

    write_csv(arrange(compareMutations, padj), here::here("TCGA_results", project, automataID, "snv_analysis",
                                              paste0("mutationProportionComparison_", project, "_", process, "_", automataID, ".csv")))


    # Pathway analysis --
    geneSelected <- mafSelected@gene.summary %>%
      dplyr::select(Hugo_Symbol, MutatedSamples)

    geneContext <- mafContext@gene.summary %>%
      dplyr::select(Hugo_Symbol, MutatedSamples)

    # Replicates the values to create a vector with all the genes repeated the right
    # number of times
    geneListSelected <- vector()
    for(i in 1:nrow(geneSelected)){
      repGeneVector <- rep(unlist(geneSelected[i, 1], use.names=FALSE), geneSelected[i,2])
      geneListSelected <- c(geneListSelected, repGeneVector)
    }

    geneListContext <- vector()
    for(i in 1:nrow(geneContext)){
      repGeneVector <- rep(unlist(geneContext[i, 1], use.names=FALSE), geneContext[i,2])
      geneListContext <- c(geneListContext, repGeneVector)
    }

    runPathwayAnalysis(project, automataID, process, geneNames = geneListSelected,
                       contextNames = geneListContext, positive.l2FC = TRUE,
                       type = "snv")



  }

  message("\n\n")
  message("=======================================================")
  message("\n\n")


}
