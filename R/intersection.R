
#' Does the intersection
#'
#' Checks which genes are differentially expressed, differentially methylated, have mutations
#' and copy number variations.
#' @importFrom readr read_csv
#' @importFrom dplyr filter pull
#' @importFrom Vennerable Venn
#' @import magrittr
#' @export


createIntersection <- function(project, automataID){

  createTCGAfolder("TCGA-OV", "eridanus")

  # Imports the gamma values to know how many processes are
  gammaValues <- read_csv(here::here("TCGA_results", project, automataID, "lpd_parameters",
                                     paste0("Sample_assigned_", project, "_", automataID, ".csv")))

  # Does a loop for each one
  processes <- unique(gammaValues$Max_LPD)

  for(process in processes){

    # Reduces the gamma values to only the selected samples
    sampleList <- gammaValues %>%
      dplyr::filter(Max_LPD == process)

    # Imports the DE genes
    DEgenes <- read_csv(here::here("TCGA_results", project, automataID, "expression_analysis",
                                   paste0("DE-Genes_", project, "_", process, "_", automataID, ".csv"))) %>%
      dplyr::filter(abs(log2FoldChange) > 1, padj <= 0.05)

    DEgenes.over <- DEgenes %>%
      filter(log2FoldChange > 0) %>%
      pull(Gene)

    DEgenes.under <- DEgenes %>%
      filter(log2FoldChange < 0) %>%
      pull(Gene)

    # Imports the DM genes
    DMgenes <- read_csv(here::here("TCGA_results", project, automataID, "methylation_analysis",
                                   paste0("DM_signifGenes_", project, "_", process, "_", automataID, ".csv")))
    if(nrow(DMgenes) > 1){

      DMgenes <- DMgenes %>%
        dplyr::filter(abs(meanlogFC) > 1)

      if(nrow(DMgenes) > 1){

        DMgenes.hyper <- DMgenes %>%
          dplyr::filter(meanlogFC > 0) %>%
          pull(Gene_Symbol)

        DMgenes.hypo <- DMgenes %>%
          dplyr::filter(meanlogFC < 0) %>%
          pull(Gene_Symbol)

      }

    } else{

      DMgenes.hyper <- DMgenes
      DMgenes.hypo <- DMgenes
    }


    # Imports the CNV genes
    aberratedGenes <- read_csv(here::here("TCGA_results", project, automataID, "cnv_analysis",
                                   paste0("aberrationsGene_", project, "_", process, "_", automataID, ".csv")))

    aberratedGenes.amp <- aberratedGenes %>%
      dplyr::filter(Amp > 0) %>%
      pull(hgnc_symbol)

    aberratedGenes.del <- aberratedGenes %>%
      dplyr::filter(Del > 0) %>%
      pull(hgnc_symbol)

    # Imports the SNV genes
    mutatedGenes <- read_csv(here::here("TCGA_results", project, automataID, "snv_analysis",
                                          paste0("mutatedGenes_", project, "_", process, "_", automataID, ".csv"))) %>%
      pull(Hugo_Symbol)



    # Creates venn diagram
    venn.over <- Venn(list(DEgenes = DEgenes.over,
                           Hypomethylated = DMgenes.hypo,
                           Amplified = aberratedGenes.amp))

    venn.under <- Venn(list(DEgenes = DEgenes.under,
                            Hypermethylated = DMgenes.hyper,
                            Deleted = aberratedGenes.del,
                            Mutated = mutatedGenes))


    pdf(here::here("TCGA_results", project, automataID, "intersect",
                   paste0("vennPlot_Overexpressed_", process, "_", automataID, ".pdf")),
        width = 10, height = 10)
    plot(venn.over, doWeights = FALSE)
    dev.off()

    save(venn.over, file = paste0("TCGA_results/", project, "/", automataID,
                                             "/report_plots/venn.over_plot_",
                                  project, "_",  process, "_", automataID, ".Rdata"))

    pdf(here::here("TCGA_results", project, automataID, "intersect",
                   paste0("vennPlot_Underexpressed_", process, "_", automataID, ".pdf")),
        width = 10, height = 10)
    plot(venn.under, type = "ellipses")
    dev.off()

    save(venn.under, file = paste0("TCGA_results/", project, "/", automataID,
                                  "/report_plots/venn.under_plot_",
                                  project, "_",  process, "_", automataID, ".Rdata"))

    # Creates summary dataset ----------------

    # Each vector can have different size so I need to fill with NAs
    list.over <- venn.over@IntersectionSets
    max.length <- max(sapply(list.over, length))
    list.over <- lapply(list.over, function(x) {c(x, rep(NA, max.length - length(x)))})

    # Now I can merge it into a dataframe
    names(list.over) <- c("NoOverlap", "DEG", "Hypomethylated", "DEG+Hypomethylated",
                          "Amplified", "DEG+Amplified", "Hypomethylated+Amplified",
                          "DEG+Hypomethylated+Amplified")

    overexpressedGenesSummary <- as.data.frame(do.call(cbind, list.over))

    #-------
    list.under <- venn.under@IntersectionSets
    max.length <- max(sapply(list.under, length))
    list.under <- lapply(list.under, function(x) {c(x, rep(NA, max.length - length(x)))})

    names(list.under) <- c("NoOverlap", "DEG", "Hypermethylated", "DEG+Hypermethylated",
                           "Deleted", "DEG+Deleted", "Hypermethylated+Deleted",
                           "DEG+Hypermethylated+Deleted", "Mutated", "DEG+Mutated",
                           "Hypermethylated+Mutated", "DEG+Hypermethylated+Mutated",
                           "Deleted+Mutated", "DEG+Deleted+Mutated",
                           "Hypermethylated+Deleted+Mutated",
                           "DEG+Hypermethylated+Deleted+Mutated")

    underexpressedGenesSummary <- as.data.frame(do.call(cbind, list.under))
    #---------

    write.csv(overexpressedGenesSummary,
              here::here("TCGA_results", project, automataID, "intersect",
                                                    paste0("overexpressedGenesSummary_", process, "_", automataID, ".csv")))

    write.csv(underexpressedGenesSummary,
              here::here("TCGA_results", project, automataID, "intersect",
                         paste0("underexpressedGenesSummary_", process, "_", automataID, ".csv")))


  }



}
