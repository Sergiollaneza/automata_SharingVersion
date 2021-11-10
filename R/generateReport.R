#' Generates report to check results
#'
#' Creates and HTML file containing the resuls from the whole automata pipeline
#' @param project Name of a TCGA project (eg. TCGA-PRAD or TCGA-BRCA).
#' @param automataID ID assigned to the analysis so the files from different analysis are not mixed.
#' @importFrom readr read_csv
#' @importFrom rmarkdown render
#' @importFrom dplyr mutate_at
#' @import here
#' @import flexdashboard
#' @export
generateReport <- function(project, automataID, microarray, workingPath = getwd()){


message("\n------GENERATING REPORT-------\n")
# Preparations ------------------------------------------------------------

  # Assigns the right method
  if(microarray == TRUE) method = "microarray" else method = "RNA-seq"

  # Creates required folders
  message("\no Creating all necessary folders...")
  createTCGAfolder(project, automataID)

  # Creates RMD file for the Main Document
  message("\no Creating RMD files for the main document...")
  createMainRMD(project, automataID, workingPath)
  # Creates RMD file for the overview valueboxes
  message("\no Creating RMD files for the valueboxes...")
  createOverviewRMD(project, automataID, workingPath)

# 01. Barcodes ------------------------------------------------------------
  message("\no Generating main document variables...")
  message("\noo Barcodes page")
  show.barcodes <- FALSE
  barcodes.data <- NULL

  # Check if the barcodes data were generated and imports them
  if(file.exists(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/TCGA_data/",
                        "Barcodes_", project, "_", automataID, ".csv"))){

    show.barcodes <- TRUE

    barcodes.data <- read_csv(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/TCGA_data/",
                                     "Barcodes_", project, "_", automataID, ".csv"))
  }


# 02. Top genes ------------------------------------------------------------
  message("\noo Genes page")
  show.genes <- FALSE
  genes.data <- NULL

  if(file.exists(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/processed_expression/",
                        "topGenes_", project, "_", automataID, ".csv"))){

    show.genes <- TRUE

    genes.data <- colnames(read_csv(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/processed_expression/",
                                     "topGenes_", project, "_", automataID, ".csv"))[,-1])


    divisor <- ceiling(length(genes.data)/5)

    genes.dataTable <- data.frame("First" = genes.data[1:divisor],
                                  "Second" = genes.data[((divisor) + 1):(divisor*2)],
                                  "Third" = genes.data[((divisor * 2) + 1):(divisor * 3)],
                                  "Fourth" = genes.data[((divisor * 3) + 1):(divisor * 4)],
                                  "Fifth" = genes.data[((divisor * 4) + 1):(divisor * 5)])




  }

# 03. LPD parameters ------------------------------------------------------------
  message("\noo LPD process page")
  show.LPD <- FALSE
  processesAndSigmas.data <- NULL
  likelihood.plot <- NULL
  bestCombination.data <- NULL
  correlations.data <- NULL
  sampleAssigned.data <- NULL
  gammaSample.plot <- NULL


  if(file.exists(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/lpd_parameters/",
                        "bestProcessesAndSigmas_", project, "_", automataID, ".csv")) &&
     file.exists(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/lpd_parameters/",
                        "bestCombination_", project, "_", automataID, ".csv")) &&
     file.exists(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/lpd_parameters/",
                        "bestRuns_", project, "_", automataID, ".csv")) &&
     file.exists(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/lpd_parameters/",
                        "correlations_", project, "_", automataID, ".csv")) &&
     file.exists(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/lpd_parameters/",
                        "Sample_assigned_", project, "_", automataID, ".csv")) &&
     file.exists(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/report_plots/",
                        "Likelihood_Plot_", project, "_", automataID, ".Rdata")) &&
     file.exists(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/report_plots/",
                        "gammaSample_Plot_", project, "_", automataID, ".Rdata"))){

    # Gene data is available
    show.LPD <- TRUE

    # Read CSV to show in the report

    # Selected processes and sigma
    processAndSigmas.data  <- read_csv(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/lpd_parameters/",
                                                "bestProcessesAndSigmas_", project, "_", automataID, ".csv")) %>%
      mutate_at(4:7, round, 4)

    # Best combination selected
    bestCombination.data  <- read_csv(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/lpd_parameters/",
                                             "bestCombination_", project, "_", automataID, ".csv"))
    # Pearson results
    correlations.data  <- read_csv(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/lpd_parameters/",
                                          "correlations_", project, "_", automataID, ".csv")) %>%
      mutate_at(3:6, round, 4) %>%
      dplyr::select(-Folder) %>%
      dplyr::select(n_process, sigma, everything()) %>%
      mutate(Significant = ifelse(pVal <= 0.05, "*", ""))

    # Gamma assigned samples
    sampleAssigned.data  <- read_csv(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/lpd_parameters/",
                                            "Sample_assigned_", project, "_", automataID, ".csv"))


    # Loads plot of Likelihood
    load(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/report_plots/",
                "Likelihood_Plot_", project, "_", automataID, ".Rdata"))

    # Loads plot of gamma values
    load(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/report_plots/",
                "gammaSample_Plot_", project, "_", automataID, ".Rdata"))


  }

  # Creates RMD file for the child pages
  message("\no Creating RMD files for thechild documents...")

  # Creates a vector that has one true for each existing LPD group
  subtypesVector <- rep(TRUE, bestCombination.data$n_process)
  subtypesVector <- c(subtypesVector, rep(FALSE, 15- bestCombination.data$n_process))

  createChildsRMD(project, automataID, workingPath, LPDvector = subtypesVector)


# 04. Signature analysis ------------------------------------------------------------
  message("\noo Checking how many processes were found...")
  KMplot.List <- list()
  DEgenes.List <- list()
  DMgenes.List <- list()
  DMsamplePlot.list <- list()
  DMchromosomePlot.list <- list()
  mutatedGenes.List <- list()
  oncoplot.List <- list()
  snpTypePlot.List <- list()
  variantClassificationPlot.List <- list()
  variantSample.List <- list()
  variantType.List <- list()
  aberratedGenesAmp.List <- list()
  aberratedGenesDel.List <- list()
  aberrationSample.List <- list()
  aberrationChromosomeAmp.List <- list()
  aberrationChromosomeDel.List <- list()
  intersectOverTable.List <- list()
  intersectUnderTable.List <- list()
  vennOver.List <- list()
  vennUnder.List <- list()
  volcanoPlot.List <- list()


  message("\nooo Creating pages for each process...")
  for(process in 1:length(subtypesVector)){

    if(subtypesVector[process] == TRUE){


      # A. Clinical plots ###
      if(file.exists(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/report_plots/",
                            "survivalOverall_plot_", project, "_LPD_", process, "_", automataID, ".Rdata"))){

        KMplot.List[[process]] <- get(load(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/report_plots/",
                                                  "survivalOverall_plot_", project, "_LPD_", process, "_", automataID, ".Rdata")))


      }

      # B. Differential expressed genes ###
      if(file.exists(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/expression_analysis/",
                            "DE-Genes_", project, "_LPD_", process, "_", automataID, ".csv"))){

        DEgenes.List[[process]] <- read_csv(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/expression_analysis/",
                                                   "DE-Genes_", project, "_LPD_", process, "_", automataID, ".csv")) %>%
          dplyr::select(Gene, log2FoldChange, padj) %>%
          mutate(Significant = ifelse(padj <= 0.05, "*", "")) %>%
          mutate_at(2:3, round, 4) %>%
          dplyr::filter(!is.na(log2FoldChange),
                        !is.na(padj))

      }

      # C. Differential methylated genes ###
      if(file.exists(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/methylation_analysis/",
                            "DM_signifGenes_", project, "_LPD_", process, "_", automataID, ".csv")) &&
         file.exists(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/report_plots/",
                            "DmChromosome_Plot_", project, "_LPD_", process, "_", automataID, ".Rdata")) &&
         file.exists(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/report_plots/",
                            "dmSample_Plot_", project, "_LPD_", process, "_", automataID, ".Rdata"))){

        DMgenes.List[[process]] <- read_csv(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/methylation_analysis/",
                                                   "DM_signifGenes_", project, "_LPD_", process, "_", automataID, ".csv"))
        if(nrow(DMgenes.List[[process]]) > 0){

          DMgenes.List[[process]] <- DMgenes.List[[process]] %>%
            mutate(meanlogFC = round(meanlogFC, 4)) %>%
            dplyr::filter(!is.na(meanlogFC)) %>%
            arrange(Chromosome)

        }


        DMsamplePlot.list[[process]] <- get(load(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/report_plots/",
                                             "dmSample_Plot_", project, "_LPD_", process, "_", automataID, ".Rdata")))

        DMchromosomePlot.list[[process]] <- get(load(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/report_plots/",
                                                            "DmChromosome_Plot_", project, "_LPD_", process, "_", automataID, ".Rdata")))

      }


      # D. SNVs ##
      if(file.exists(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/snv_analysis/",
                            "mutationProportionComparison_", project, "_LPD_", process, "_", automataID, ".csv")) &&
         file.exists(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/report_plots/",
                            "oncoplot_Plot_", project, "_LPD_", process, "_", automataID, ".Rdata")) &&
         file.exists(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/report_plots/",
                            "snpType_Plot_", project, "_LPD_", process, "_", automataID, ".Rdata")) &&
         file.exists(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/report_plots/",
                            "variantClassification_Plot_", project, "_LPD_", process, "_", automataID, ".Rdata")) &&
         file.exists(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/report_plots/",
                            "variantSample_Plot_", project, "_LPD_", process, "_", automataID, ".Rdata")) &&
         file.exists(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/report_plots/",
                            "variantType_Plot_", project, "_LPD_", process, "_", automataID, ".Rdata"))){

        mutatedGenes.List[[process]] <- read_csv(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/snv_analysis/",
                                                        "mutationProportionComparison_", project, "_LPD_", process, "_", automataID, ".csv")) %>%
          dplyr::select(-c(NotMutated.Selected, NotMutated.Context)) %>%
          mutate(Significant = ifelse(padj <= 0.05, "*", ""))


        oncoplot.List[[process]] <- get(load(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/report_plots/",
                                                        "oncoplot_Plot_", project, "_LPD_", process, "_", automataID, ".Rdata")))

        snpTypePlot.List[[process]] <- get(load(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/report_plots/",
                                                            "snpType_Plot_", project, "_LPD_", process, "_", automataID, ".Rdata")))

        variantClassificationPlot.List[[process]] <- get(load(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/report_plots/",
                                                    "variantClassification_Plot_", project, "_LPD_", process, "_", automataID, ".Rdata")))

        variantSample.List[[process]] <- get(load(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/report_plots/",
                                                       "variantSample_Plot_", project, "_LPD_", process, "_", automataID, ".Rdata")))

        variantType.List[[process]] <- get(load(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/report_plots/",
                                                         "variantType_Plot_", project, "_LPD_", process, "_", automataID, ".Rdata")))


      }

      # CNVs
      if(file.exists(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/cnv_analysis/",
                            "aberrationsGene_", project, "_LPD_", process, "_", automataID, ".csv")) &&
         file.exists(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/report_plots/",
                            "aberrationSample_plot_", project, "_LPD_", process, "_", automataID, ".Rdata")) &&
         file.exists(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/report_plots/",
                             "aberrationChromosome_Amp_plot_", project, "_LPD_", process, "_", automataID, ".Rdata")) &&
         file.exists(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/report_plots/",
                            "aberrationChromosome_Del_plot_", project, "_LPD_", process, "_", automataID, ".Rdata"))){

        aberratedGenesAmp.List[[process]] <- read_csv(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/cnv_analysis/",
                                                          "aberrationsGene_", project, "_LPD_", process, "_", automataID, ".csv")) %>%
          dplyr::select(-c(notAmp.Selected, notAmp.Context, Del.Selected, Del.Context, notDel.Selected, notDel.Context, Del.chisq, Del.padj)) %>%
          filter(hgnc_symbol != "Unknown")

        aberratedGenesDel.List[[process]] <- read_csv(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/cnv_analysis/",
                                                             "aberrationsGene_", project, "_LPD_", process, "_", automataID, ".csv")) %>%
          dplyr::select(-c(notDel.Selected, notDel.Context, Amp.Selected, Amp.Context, notAmp.Selected, notAmp.Context, Amp.chisq, Amp.padj)) %>%
          filter(hgnc_symbol != "Unknown")

        aberrationSample.List[[process]] <- get(load(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/report_plots/",
                                                            "aberrationSample_plot_", project, "_LPD_", process, "_", automataID, ".Rdata")))

        aberrationChromosomeAmp.List[[process]] <- get(load(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/report_plots/",
                                                                   "aberrationChromosome_Amp_plot_", project, "_LPD_", process, "_", automataID, ".Rdata")))

        aberrationChromosomeDel.List[[process]] <- get(load(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/report_plots/",
                                                                "aberrationChromosome_Del_plot_", project, "_LPD_", process, "_", automataID, ".Rdata")))


      }

      # Intersect
      if(file.exists(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/intersect/",
                            "overexpressedGenesSummary", "_LPD_", process, "_", automataID, ".csv")) &&
         file.exists(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/intersect/",
                            "underexpressedGenesSummary", "_LPD_", process, "_", automataID, ".csv")) &&
         file.exists(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/report_plots/",
                            "venn.over_plot_", project, "_LPD_", process, "_", automataID, ".Rdata"))){

        intersectOverTable.List[[process]] <- read_csv(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/intersect/",
                                                              "overexpressedGenesSummary", "_LPD_", process, "_", automataID, ".csv")) %>%
          dplyr::select(-c(X1, NoOverlap, DEG, Hypomethylated, Amplified))

        intersectUnderTable.List[[process]] <- read_csv(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/intersect/",
                                                                    "underexpressedGenesSummary", "_LPD_", process, "_", automataID, ".csv")) %>%
          dplyr::select(-c(X1, NoOverlap, DEG, Hypermethylated, Deleted, Mutated))

        colnames(intersectUnderTable.List[[process]]) <- c("EH", "EA", "HA", "EHA", "EM", "HM", "EHM", "AM", "EAM", "HAM", "EHAM")

        vennOver.List[[process]] <- get(load(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/report_plots/",
                                                    "venn.over_plot_", project, "_LPD_", process, "_", automataID, ".Rdata")))
        vennUnder.List[[process]] <- get(load(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/report_plots/",
                                                     "venn.under_plot_", project, "_LPD_", process, "_", automataID, ".Rdata")))

      }


      generateSubChilds(project, automataID, process, workingPath)

    }

  }

  # Some objects are overall and don't need the loop ----

  # Clinical plots
  if(file.exists(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/report_plots/",
                        "survivalOverall_plot_", project, "_", automataID, ".Rdata"))){

    KMplot.List[[16]] <- get(load(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/report_plots/",
                                         "survivalOverall_plot_", project, "_", automataID, ".Rdata")))


  }

  # Differential expressed genes
  if(file.exists(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/report_plots/",
                        "deGenes_Plot_", project, "_", automataID, ".Rdata"))){


    volcanoPlot.List <- get(load(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/report_plots/",
                                        "deGenes_Plot_", project, "_", automataID, ".Rdata")))

  }


  # 05. Overall cancer analysis ------------------------------------------------------------
  batchEffect.data <- NULL
  raceProportion.data <- NULL


  if(file.exists(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/batch_effect/",
                          "sampleTypeProportions_", project, "_", automataID, ".csv"))){

      batchEffect.data <- read_csv(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/batch_effect/",
                                         "sampleTypeProportions_", project, "_", automataID, ".csv")) %>%
        mutate_at(vars(starts_with("chisq")), round, 4)


 }

  if(file.exists(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/clinical_analysis/",
                        "raceProportions_", project, "_", automataID, ".csv"))){

    raceProportion.data <- read_csv(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/clinical_analysis/",
                                        "raceProportions_", project, "_", automataID, ".csv")) %>%
      mutate_at(vars(starts_with("chisq")), round, 4)


  }


  if(file.exists(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/dendrogram/",
                          "Something_Plot_", project, "_LPD_", process, "_", automataID, ".Rdata"))){

      print("Obviously, not ready")


  }


  # 06. Pathway analysis ------------------------------------------------------------
  message("\noo Checking how many processes were found...")
  DEgoOver.List <- list()
  DEgoUnder.List <- list()
  DEkeggOver.List <- list()
  DEkeggUnder.List <- list()
  DEHallmarkOver.List <- list()
  DEHallmarkUnder.List <- list()



  message("\nooo Creating pages for each process...")
  for(process in 1:length(subtypesVector)){

    if(subtypesVector[process] == TRUE){


      if(file.exists(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/expression_analysis/pathway_analysis/",
                            "egoAnalysis_overexpressedGenes_LPD_", process, "_", automataID, ".csv"))){

        DEgoOver.List[[process]] <- read_csv(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/expression_analysis/pathway_analysis/",
                                                    "egoAnalysis_overexpressedGenes_LPD_", process, "_", automataID, ".csv")) %>%
          dplyr::select(-c(BgRatio, pvalue, qvalue, geneID, Count)) %>%
          dplyr::filter(p.adjust <= 0.05)

        DEgoUnder.List[[process]] <- read_csv(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/expression_analysis/pathway_analysis/",
                                                     "egoAnalysis_underexpressedGenes_LPD_", process, "_", automataID, ".csv")) %>%
          dplyr::select(-c(BgRatio, pvalue, qvalue, geneID, Count)) %>%
          dplyr::filter(p.adjust <= 0.05)

        DEkeggOver.List[[process]] <- read_csv(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/expression_analysis/pathway_analysis/",
                                                      "keggAnalysis_overexpressedGenes_LPD_", process, "_", automataID, ".csv")) %>%
          dplyr::select(-c(BgRatio, pvalue, qvalue, geneID, Count)) %>%
          dplyr::filter(p.adjust <= 0.05)

        DEkeggUnder.List[[process]] <- read_csv(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/expression_analysis/pathway_analysis/",
                                                       "keggAnalysis_underexpressedGenes_LPD_", process, "_", automataID, ".csv")) %>%
          dplyr::select(-c(BgRatio, pvalue, qvalue, geneID, Count)) %>%
          dplyr::filter(p.adjust <= 0.05)

        DEHallmarkOver.List[[process]] <- read_csv(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/expression_analysis/pathway_analysis/",
                                                          "hallmarkAnalysis_overexpressedGenes_LPD_", process, "_", automataID, ".csv")) %>%
          dplyr::select(-c(BgRatio, pvalue, qvalue, Count, entrezID, hgnc_symbol)) %>%
          dplyr::filter(p.adjust <= 0.05)

        DEHallmarkUnder.List[[process]] <- read_csv(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/expression_analysis/pathway_analysis/",
                                                           "hallmarkAnalysis_underexpressedGenes_LPD_", process, "_", automataID, ".csv")) %>%
          dplyr::select(-c(BgRatio, pvalue, qvalue, Count, entrezID, hgnc_symbol)) %>%
          dplyr::filter(p.adjust <= 0.05)


      }

      generatePathways(project, automataID, process, workingPath)

    }

  }


# Render main document ----------------------------------------------------
  message("\no Render document")
  rmarkdown::render(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/Report_RMD/mainDocument.Rmd"),
                    output_file = paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/Report/Automata_Report_01.html"),
                    params = list(
                      project = project,
                      automataID = automataID,
                      method = method,
                      path = workingPath,
                      show_barcodes = show.barcodes,
                      barcodes = barcodes.data,
                      show_genes = show.genes,
                      genes = genes.dataTable,
                      show_LPD = show.LPD,
                      processAndSigmas = processAndSigmas.data,
                      likelihood_plot = likelihood.plot,
                      bestCombination = bestCombination.data,
                      correlations = correlations.data,
                      sampleAssigned = sampleAssigned.data,
                      gammaSample_plot = gammaSample.plot,
                      KMplot_List = KMplot.List,
                      DEgenes_List = DEgenes.List,
                      volcanoPlot_List = volcanoPlot.List,
                      DMgenes_List = DMgenes.List,
                      DMsamplePlot_List = DMsamplePlot.list,
                      DMchromosomePlot_List = DMchromosomePlot.list,
                      mutatedGenes_List = mutatedGenes.List,
                      oncoplot_List = oncoplot.List,
                      snpTypePlot_List = snpTypePlot.List,
                      variantClassificationPlot_List = variantClassificationPlot.List,
                      variantSample_List = variantSample.List,
                      variantType_List = variantType.List,
                      aberrationSample_List = aberrationSample.List,
                      intersectOverTable_List = intersectOverTable.List,
                      intersectUnderTable_List = intersectUnderTable.List,
                      batchEffect = batchEffect.data,
                      raceProportion = raceProportion.data,
                      DEgoOver_List = DEgoOver.List,
                      DEgoUnder_List = DEgoUnder.List,
                      DEkeggOver_List = DEkeggOver.List,
                      DEkeggUnder_List = DEkeggUnder.List,
                      DEHallmarkOver_List = DEHallmarkOver.List,
                      DEHallmarkUnder_List = DEHallmarkUnder.List,
                      aberratedGenesAmp_List = aberratedGenesAmp.List,
                      aberratedGenesDel_List = aberratedGenesDel.List,
                      aberrationChromosomeAmp_List = aberrationChromosomeAmp.List,
                      aberrationChromosomeDel_List = aberrationChromosomeDel.List,
                      vennUnder_List = vennUnder.List,
                      vennOver_List = vennOver.List))




}
