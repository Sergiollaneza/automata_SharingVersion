#' Analysis of methylation data
#'
#' Performs an analysis of the methylation data downloaded from the TCGA comparing each LPD group.
#' @param project Name of the project to analyse according to the TCGA (e.g. TCGA-PRAD)
#' @param automataID Custom ID tag
#' @import magrittr
#' @import here
#' @import ggplot2
#' @importFrom readr read_csv write_csv
#' @importFrom dplyr filter arrange mutate distinct summarise group_by
#' @importFrom Biobase AnnotatedDataFrame ExpressionSet
#' @importFrom limma lmFit eBayes topTable
#' @importFrom gridExtra tableGrob grid.arrange
#' @importFrom tidyr separate_rows
#' @importFrom tibble column_to_rownames
#' @export

methylationAnalysis <- function(project, automataID){

  message("\n\n o METHYLATION ANALYSIS ---------------------\n\n")

  # Import data -------------------------------------------------------------
  message("oo Importing the data....")
  # Imports to which LPD group is each sample assigned to
  gammaValues <- read_csv(here::here("TCGA_results", project, automataID, "lpd_parameters",
                                     paste0("Sample_assigned_", project, "_", automataID, ".csv"))) %>%
    mutate(Sample = substr(Aliquot, 1, 16))

  # Extracts all the processes found for the dataset
  processes <- unique(gammaValues$Max_LPD)

  # Imports the beta values downloaded from the TCGA and change the sample code to the standard one
  betaValues <- read_csv(here::here("TCGA_results", project, automataID, "TCGA_data",
                                    paste0("meth_matrix_", project, "_", automataID, ".csv")))


  message("ooo Cleaning the data....")
  # Cleans duplicated samples
  sampleNames <- data.frame("Barcode" = names(betaValues)[11:ncol(betaValues)]) %>%
    separate(Barcode, c("Project", "TSS", "Participant", "Sample_Vial", "Portion_Analyte", "Plate", "Center"), "-", remove = FALSE) %>%
    separate(Sample_Vial, c("Sample", "Vial"), "(?<=[0-9])(?=[A-Z])") %>%
    separate(Portion_Analyte, c("Portion", "Analyte"), "(?<=[0-9])(?=[A-Z])") %>%
    group_by(TSS, Participant, Sample) %>%
    arrange(Vial, desc(Plate)) %>%
    slice(1) %>%
    ungroup() %>%
    select(-c(Project, TSS, Participant, Sample, Vial, Portion, Analyte, Plate, Center)) %>%
    pull()

  beta.matrix <- betaValues[,11:ncol(betaValues)][colnames(betaValues[,11:ncol(betaValues)]) %in% sampleNames]
  names(beta.matrix) <- gsub("\\.", "-", substr(names(beta.matrix), 1, 16))

  # Calculates differential methylated genes -------------------------------
  message("oo Calculating differentially methylated genes....")
  # To use limma to calculate it, first we need to create an expressionset

  # Creates the assay and feature data needed to run limma
  message("ooo Creating assay and feature object")
  assay <- as.matrix(beta.matrix)
  assay <- log2(assay/(1 - assay))

  feature <- AnnotatedDataFrame(betaValues[,1:10])

  # Generates methylation level per gene file -------------------------------
  p <- 0
  probeSeparated.list <- list()

  # Loops for each process changing the selected process vs all others
  for(process in processes){
    message("ooo Calculating methylation level for ", process)

    p <- p + 1

    # Assigns a binary (1/0) if the sample belongs to the right process
    selectedSamples <- gammaValues %>%
      separate(Aliquot, c("Project", "TSS", "Participant", "Sample_Vial", "Portion_Analyte", "Plate", "Center"), "-", remove = FALSE) %>%
      separate(Sample_Vial, c("Sample", "Vial"), "(?<=[0-9])(?=[A-Z])") %>%
      separate(Portion_Analyte, c("Portion", "Analyte"), "(?<=[0-9])(?=[A-Z])") %>%
      group_by(TSS, Participant, Sample) %>%
      arrange(Vial, desc(Plate)) %>%
      dplyr::slice(1) %>%
      ungroup() %>%
      select(-c(Project, TSS, Participant, Sample, Vial, Portion, Analyte, Plate, Center)) %>%
      mutate(Sample = substr(Aliquot, 1, 16)) %>%
      filter(Sample %in% substr(sampleNames, 0, 16)) %>%
      dplyr::mutate(Selected = as.factor(ifelse(Max_LPD == process, 1, 0))) %>%
      mutate(Order = match(Sample, colnames(assay))) %>%
      arrange(Order) %>%
      select(-c(Order)) %>%
      select(-Aliquot) %>%
      as.data.frame()

    message("ooo Creating pehnotype object for ", process)
    # Builds the phenotype data needed to run limma
    rownames(selectedSamples) <- selectedSamples$Sample
    pheno <- AnnotatedDataFrame(selectedSamples)

    # Builds the set with the assay, pehnotype and feature data
    set <- ExpressionSet(assayData = assay,
                         phenoData = pheno,
                         featureData = feature)

    if(length(levels(pData(set)$Selected)) < 2){
      probeSeparated.list[[p]] <- NA
      next()
    }

    designMatrix <- model.matrix(~Selected, data = pData(set))



    # Calculates probe region using limma
    message("ooo Fitting the model for ", process)
    fit <- lmFit(set, designMatrix)
    fit2 <- eBayes(fit)
    probe <- topTable(fit2, adjust="BH",  num=Inf)

    # Tidies the output
    probeSeparated <- probe %>%
      separate_rows(Gene_Symbol, Gene_Type, Transcript_ID, Position_to_TSS, sep = ";") %>%
      distinct() %>%
      filter(Gene_Symbol != ".") %>%
      dplyr::rename(Composite.Element.REF = X1) %>%
      arrange(Composite.Element.REF) %>%
      filter(Chromosome != "*") %>%
      mutate(Chromosome = substr(Chromosome, 4, 5),
             Chromosome = gsub("X", 23, Chromosome),
             Chromosome = gsub("Y", 24, Chromosome))

    probeSeparated.list[[p]] <- probeSeparated

    write_csv(probeSeparated, here::here("TCGA_results", project, automataID, "methylation_analysis",
                                         paste0("Probe_Region_", project, "_", process, "_", automataID, ".csv")))

  }

  # Calculates significant differentially methylated genes --------------------------

  p <- 0

  for (process in processes){

    p <- p + 1

    if(is.na(probeSeparated.list[[p]])){
      next()
    }

    # Keeps only each gene once taking the average log2FC and only the statistical significant ones
    message("ooo Summarising multiple matches for each gene for ", process)
    uniqueGenes <- probeSeparated.list[[p]] %>%
      filter(adj.P.Val <= 0.05) %>%
      distinct(Chromosome, Gene_Symbol, Gene_Type, logFC, .keep_all = TRUE) %>%
      group_by(Chromosome, Gene_Symbol) %>%
      summarise(meanlogFC = mean(logFC, na.rm = TRUE), p.val = dplyr::first(adj.P.Val)) %>%
      dplyr::filter(abs(meanlogFC) >= 1.5) %>%
      arrange(desc(meanlogFC))

    write_csv(uniqueGenes, here::here("TCGA_results", project, automataID, "methylation_analysis",
                                      paste0("DM_SignifGenes_", project, "_", process, "_", automataID, ".csv")))


    if(nrow(uniqueGenes) > 0){
      message("oo Runs pathway analysis for ", process)
      # Performs a pathway analysis
      hypermeth.genes <- pull(filter(uniqueGenes, meanlogFC >= 1), Gene_Symbol)

      if(length(hypermeth.genes) > 0){

        runPathwayAnalysis(project, automataID, process, geneNames = hypermeth.genes,
                           contextNames = unique(probeSeparated.list[[p]]$Gene_Symbol),
                           positive.l2FC = TRUE, type = "methylation")
      }

      hypometh.genes <- pull(filter(uniqueGenes, meanlogFC <= -1), Gene_Symbol)

      if(length(hypometh.genes) > 0){

        runPathwayAnalysis(project, automataID, process, geneNames = hypometh.genes,
                           contextNames = unique(probeSeparated.list[[p]]$Gene_Symbol),
                           positive.l2FC = FALSE, type = "methylation")

      }

    }



  }

  # Calculate ratio hyper/hypo per sample --------------------------
  message("oo Calculating the ratio hypermethylation/hypomethylation for each sample")
  # Creates a list object with the number of 0 and 1 for each sample, then is converted to DF
  betaSample <- betaValues[11:ncol(betaValues)] %>%
    mutate_all(~ case_when(. >= 0.8 ~ 1,
                           . <= 0.2 ~ 0)) %>%
    purrr::map(~table(.)) %>%
    bind_cols() %>%
    t() %>%
    as.data.frame() %>%
    mutate(Sample = substr(rownames(.), 1, 16),
           Hypomethylated = V1,
           Hypermethylated = V2,
           Ratio = Hypermethylated / Hypomethylated,
           Percentage.Hypermethylated = Hypermethylated * 100 / (Hypermethylated + Hypomethylated),
           Percentage.Hypomethylated = 100 - Percentage.Hypermethylated) %>%
    inner_join(gammaValues, by = "Sample") %>%
    dplyr::select(-c(V1, V2))

  write_csv(betaSample, here::here("TCGA_results", project, automataID, "methylation_analysis",
                                   paste0("DMSample_", project, "_", automataID, ".csv")))

  betaSample.pvalue.tibble<- tibble()

  for(process in processes){
    message("ooo In ", process)
    betaSample.assigned <- dplyr::mutate(betaSample, Selected = as.factor(ifelse(Max_LPD == process, process, "Other_LPD")))

    if(length(levels(betaSample.assigned$Selected)) < 2){
      message("oooo Only one level found, skipping...")
      next()
    }

    betaSample.pvalue <- wilcox.test(formula = Ratio ~ Selected, data = betaSample.assigned)$p.value
    betaSample.pvalue.tibble <- bind_rows(betaSample.pvalue.tibble,
                                          setNames(c(process, betaSample.pvalue), c("Process", "Pvalue")))


    plot <- ggplot(betaSample.assigned,
                   aes(x = Selected, y = Ratio,  fill = Selected)) +
      geom_violin(trim = FALSE) +
      coord_flip() +
      theme_minimal() +
      geom_boxplot(width=0.1) +
      scale_fill_manual(values = c("#81ab72", '#4e7741')) +
      ggtitle("Hyper/hypo methylation ratio per sample",
              subtitle = paste0("p-value: ", round(betaSample.pvalue, digits = 5))) +
      ylim(0, 2) +
      theme(axis.title.y = element_blank(),
            plot.title = element_text(hjust = 0.5, face = "bold"),
            legend.position = "none")

    ggsave(here::here("TCGA_results", project, automataID, "methylation_analysis",
                      paste0("DMsample_Plot_", project, "_", process, "_", automataID, ".pdf")),
           plot, width = 6, height = 4, units = "in")

    dmSample_Plot <- plot

    save(dmSample_Plot, file = paste0("TCGA_results/", project, "/", automataID,
                                      "/report_plots/dmSample_Plot_", project, "_", process, "_", automataID, ".Rdata"))
  }

  # Export pvalues
  write_csv(betaSample.pvalue.tibble, here::here("TCGA_results", project, automataID, "methylation_analysis",
                                                 paste0("DMSample_pValues_", project, "_", automataID, ".csv")))

  # Calculates ratio hyper/hypo per chromosome ------------------------------
  message("oo Calculating the ratio hypermethylation/hypomethylation for each chromosome")

  betaTransformed <- betaValues[11:ncol(betaValues)] %>%
    mutate_all(~ case_when(. >= 0.8 ~ 1,
                           . <= 0.2 ~ 0))

  colnames(betaTransformed) <- substr(colnames(betaTransformed), 1, 16)

  for(process in processes){

    # Creates a vector with the samples belonging to the current process and the others
    names.selected <- colnames(betaTransformed[,which(colnames(betaTransformed) %in% filter(gammaValues, Max_LPD == process)$Sample)])
    names.context <- colnames(betaTransformed[,which(colnames(betaTransformed) %in% filter(gammaValues, Max_LPD != process)$Sample)])

    # Creates a dataframe with all the data per cgSite
    betaChromosome <- data.frame(cgSite <- betaValues$X1,
                                 Chromosome = betaValues$Chromosome,
                                 Hypermeth.selected = rowSums(betaTransformed[,names.selected], na.rm = TRUE),
                                 Hypermeth.context = rowSums(betaTransformed[,names.context], na.rm = TRUE)) %>%
      mutate(Hypometh.selected = length(names.selected) - Hypermeth.selected - rowSums(is.na(betaTransformed[,names.selected]), na.rm = TRUE),
             Hypometh.context = length(names.context) - Hypermeth.context - rowSums(is.na(betaTransformed[,names.context]), na.rm = TRUE)) %>%
      group_by(Chromosome) %>%
      summarise(Hypermeth.selected.sum = sum(Hypermeth.selected),
                Hypometh.selected.sum = sum(Hypometh.selected),
                Hypermeth.context.sum = sum(Hypermeth.context),
                Hypometh.context.sum = sum(Hypometh.context),
                Ratio.selected = Hypermeth.selected.sum/Hypometh.selected.sum,
                Ratio.context = Hypermeth.context.sum/Hypometh.context.sum) %>%
      filter(Chromosome != "*") %>%
      mutate(Chromosome = substr(Chromosome, 4, 5),
             Chromosome = gsub("X", 23, Chromosome),
             Chromosome = gsub("Y", 24, Chromosome)) %>%
      arrange(as.numeric(Chromosome))

    # Calculates differences in proportions selected/context
    betaChromosome$chisq <- apply(betaChromosome, 1, function(x){
      tbl = matrix(as.numeric(x[2:5]), ncol = 2, byrow = TRUE)
      chisq.test(tbl)$p.value
    })

    betaChromosome$padj <- p.adjust(betaChromosome$chisq, method = "BH")
    betaChromosome <- mutate(betaChromosome, Significative = ifelse(betaChromosome$padj <= 0.05, "*", ""))

    # Export it
    write_csv(betaChromosome, here::here("TCGA_results", project, automataID, "methylation_analysis",
                                         paste0("DMchromosomes_", project, "_", process, "_", automataID, ".csv")))


    # Plots it
    betaChromosome.plot <- gather(dplyr::select(betaChromosome, Chromosome, Ratio.selected, Ratio.context), "Ratio", "Value", 2:3) %>%
      mutate(Ratio = ifelse(Ratio == "Ratio.selected", process, "Other LPD"),
             Significative = c(betaChromosome$Significative, rep("", length(betaChromosome$Significative))))

    plot <- ggplot(betaChromosome.plot, aes(x = reorder(Chromosome, as.integer(Chromosome)), y = Value, fill = Ratio,
                                            width=.5, label = Significative)) +
      geom_bar(stat = "identity", position=position_dodge()) +
      theme_minimal() +
      ggtitle("Hyper/hypo methylation ratio per chromosome",
              subtitle = paste("Mean:", round(mean(betaChromosome.plot$Value), digits = 2))) +
      ylim(0, 2) +
      geom_hline(yintercept=mean(betaChromosome.plot$Value), linetype="dashed", color = "red") +
      theme(axis.title.y = element_blank(),
            axis.title.x = element_blank(),
            plot.title = element_text(hjust = 0.5, face = "bold")) +
      scale_fill_manual(values=c("#81ab72", '#4e7741'), labels = c(process, "Other LPD")) +
      geom_text(vjust = -3, size = 6)

    ggsave(here::here("TCGA_results", project, automataID, "methylation_analysis",
                      paste0("DMchromosome_Plot_", project, "_", process, "_", automataID, ".pdf")),
           plot, width = 12, height = 6, units = "in")

    dmChromosome_Plot <- plot

    save(dmChromosome_Plot, file = paste0("TCGA_results/", project, "/", automataID,
                                          "/report_plots/dmChromosome_Plot_", project, "_", process, "_", automataID, ".Rdata"))


  }


}
