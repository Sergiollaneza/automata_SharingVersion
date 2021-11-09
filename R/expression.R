#' Get the differentially expressed genes
#'
#' Creates CSV files with the log2Folf of the top 500 genes and plots it in a tornado plot for each LPD process.
#'
#' @param project .
#' @param automataID .
#' @param microarray Defaults to FALSE
#' @importFrom readr read_csv write_csv
#' @import here
#' @import magrittr
#' @importFrom tibble rownames_to_column as_tibble
#' @importFrom dplyr mutate rename select case_when filter left_join
#' @importFrom tidyr gather spread
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results
#' @importFrom limma lmFit eBayes topTable
#' @importFrom purrr map_df
#' @import ggplot2
#' @importFrom ggpubr theme_pubr
#' @export

diffExpGenes <- function(project, automataID, microarray = FALSE){

  message("\n\no ANALYSIS OF DIFFERENTIALLY EXPRESSED GENES")

  # Creates a variable to store the summary or uses one already generated
  if(file.exists(here::here("TCGA_results", project, automataID,
                            paste0("Summary_", project, "_", automataID, ".csv")))){

    summary <- read_csv(here::here("TCGA_results", project, automataID,
                                   paste0("Summary_", project, "_", automataID, ".csv"))) %>%
      mutate(across(everything(), as.character))
  } else{
    summary <- tibble()
  }



  # Read the assigned samples
  message("oo Importing data")
  message("ooo Gamma values")
  samples <- read_csv(here::here("TCGA_results", project, automataID, "lpd_parameters",
                                 paste0("Sample_assigned_", project, "_", automataID, ".csv")))

  # Assign method
  if(microarray == TRUE){
    method = "microarray"

    message("ooo Transcriptome data is ", method)
    expressionMatrix <- read.csv(here::here("TCGA_results", project, automataID, "processed_expression",
                                            sprintf("clean_matrix_%s_%s_%s.csv", method, project, automataID))) %>%
      distinct(X, .keep_all = TRUE)

    rownames(expressionMatrix) <- expressionMatrix$X
    expressionMatrix <- expressionMatrix[,-1]
    hgncMatrix <- t(expressionMatrix)

    # Some rows has the word "Signal_" in microarray so let's get rid of it
    signalRows <- grep("Signal_", rownames(hgncMatrix))
    rownames(hgncMatrix)[signalRows] <- substr(rownames(hgncMatrix[signalRows,]), 8, 35)

  } else{
    method = "RNA-seq"
    message("ooo Transcriptome data is", method)

    expressionMatrix <- read.csv(here::here("TCGA_results", project, automataID, "processed_expression",
                                            paste0("clean_matrix_", method, "_", project, "_", automataID, ".csv")),
                                 row.names = 1)
    rownames(expressionMatrix) <- sub("\\..*","", rownames(expressionMatrix))
    expressionMatrix_t <- t(expressionMatrix)

    data(MARTannotations)

    ensemblNames <- as_tibble(colnames(expressionMatrix_t)) %>% dplyr::rename(ensembl_gene_id = "value")

    geneNames <- left_join(ensemblNames, MARTannotations) %>%
      distinct(ensembl_gene_id, .keep_all = TRUE) %>%
      mutate(hgnc_symbol = ifelse(hgnc_symbol == "", ensembl_gene_id, hgnc_symbol))

    colnames(expressionMatrix_t) <- geneNames$hgnc_symbol

    expressionMatrix_t <- expressionMatrix_t[,!(is.na(colnames(expressionMatrix_t)))]
    expressionMatrix_t <- expressionMatrix_t[,!(colnames(expressionMatrix_t) %in% c(""))]

    hgncMatrix <- expressionMatrix_t

  }

  geneData <- hgncMatrix %>%
    as.data.frame %>%
    rownames_to_column(var = "Sample") %>%
    dplyr::mutate(Sample = gsub("\\.", "-", substr(Sample, 1, 16)))  %>%
    gather(variable, value, -Sample) %>%
    spread(Sample, value) %>%
    dplyr::rename(Gene = variable)

  # Here I need to loop for each subtype
  processes <- unique(samples$Max_LPD)

  DEGenes.plot <- list()
  i <- 1

  for(process in processes){
    message("oo Running a differential analysis on ", process)

    # RNA-seq and microarray have different ways to calculate the expression.
    # RNA-seq needs to go through DESeq2 while microarray needs limma
    if(microarray != TRUE){

      # Does one subtype vs all the others
      selectedSamples <- samples %>%
        dplyr::mutate(Selected = as.factor(ifelse(Max_LPD == process, 1, 0)))

      # Runs DESeq to calculate differential expression
      des <- DESeqDataSetFromMatrix(countData = geneData[,-1],
                                    colData = selectedSamples,
                                    design = ~ Selected) %>%
        DESeq()

      # Calculates logFold
      desResults <- as.data.frame(results(des, contrast = c("Selected", 1, 0))) %>%
        mutate(Gene = geneData$Gene) %>%
        dplyr::select(Gene, everything()) %>%
        arrange(padj, log2FoldChange)

    } else{

      # Need to create an expression set ----------

      # Creates assy data that is the matrix of expression
      assay <- geneData[,-1] %>%
        map_df(as.numeric) %>%
        as.matrix()

      rownames(assay) <- geneData$Gene

      # Creates the phenotype dataframe and reorders it to match assay
      selectedSamples <- samples %>%
        dplyr::mutate(Selected = as.factor(ifelse(Max_LPD == process, 1, 0)))

      selectedSamples <- selectedSamples[order(selectedSamples$Sample, colnames(assay)),]

      pheno <- selectedSamples %>%
        column_to_rownames("Sample") %>%
        AnnotatedDataFrame()


      # Generates set
      set <- ExpressionSet(assayData = assay,
                           phenoData = pheno)

      # Run limma to calculate differentially expressed genes
      designMatrix <- model.matrix(~Selected, data = pData(set))
      fit <- lmFit(set, designMatrix)
      fit2 <- eBayes(fit)

      desResults <- topTable(fit2, adjust="BH",  num=Inf) %>%
        rownames_to_column("Gene") %>%
        dplyr::rename(log2FoldChange = "logFC", padj = "adj.P.Val")

    }

    # Creates summary ##
    nOver <- desResults %>%
      dplyr::filter(log2FoldChange >= 1, padj <= 0.05) %>%
      nrow()

    summaryVector <- c(Analysis = "Differentially expressed genes",
                       Test = "Diffentially overexpressed genes",
                       LPD_Group = process,
                       Significant = ifelse(nOver > 0, TRUE, FALSE),
                       n = nOver,
                       p.value = NA,
                       residual = NA,
                       date = as.character(Sys.Date()))

    summary <- bind_rows(summary, summaryVector)


    nUnder <- desResults %>%
      dplyr::filter(log2FoldChange <= -1, padj <= 0.05) %>%
      nrow()

    summaryVector <- c(Analysis = "Differentially expressed genes",
                       Test = "Diffentially underexpressed genes",
                       LPD_Group = process,
                       Significant = ifelse(nUnder > 0, TRUE, FALSE),
                       n = nUnder,
                       p.value = NA,
                       residual = NA,
                       date = as.character(Sys.Date()))

    summary <- bind_rows(summary, summaryVector)


    write_csv(desResults, here::here("TCGA_results", project, automataID, "expression_analysis",
                                     paste0("DE-Genes_", project, "_", process, "_", automataID, ".csv")))

    # Separates significant from not significant and picks the top 20 to show in the volcano plot
    significantGenes <- desResults %>%
      arrange(desc(abs(log2FoldChange))) %>%
      mutate(Expression = case_when(
        log2FoldChange >= 1 & padj <= 0.05 ~ "Overexpressed",
        log2FoldChange <= -1 & padj <= 0.05 ~ "Underexpressed",
        log2FoldChange > -1 & log2FoldChange < 1 & padj > 0.05 ~ "Not_differential",
        TRUE ~ "Not_significant")) %>%
      mutate(Show.gene = case_when(
        (Expression == "Underexpressed" | Expression == "Overexpressed") & row_number() < 21 ~ Gene,
        TRUE ~ ""))

    runPathwayAnalysis(project = project, automataID = automataID, process = process,
                       geneNames = pull(filter(significantGenes, Expression == "Overexpressed"), Gene),
                       contextNames = pull(significantGenes, Gene),
                       positive.l2FC = TRUE,
                       type = "expression")

    runPathwayAnalysis(project = project, automataID = automataID, process = process,
                       geneNames = pull(filter(significantGenes, Expression == "Underexpressed"), Gene),
                       contextNames = pull(significantGenes, Gene),
                       positive.l2FC = FALSE,
                       type = "expression")

    # Plots a volcano plot
    p <- ggplot(significantGenes, aes(x = log2FoldChange, y = -log10(padj), text = paste("Gene:", Gene), label = Show.gene)) +
      geom_point(aes(colour = Expression), size = 1.6) +
      scale_color_manual(values = c("Overexpressed"= "#0083ff", "Underexpressed"="#173bd8",
                                    "Not_significant"= "#576293", "Not_differential"= "#3d4568")) +
      geom_vline(xintercept = c(1, -1), linetype="dashed") +
      geom_hline(yintercept = -log10(0.05), linetype="dashed") +
      theme_bw() +
      theme(legend.position = "top") +
      geom_text(hjust=0, vjust=0, size = 2.5) +
      ggtitle(paste("Differentially expressed genes for", process, project, automataID, sep = " "))

    ggsave(here::here("TCGA_results", project, automataID, "expression_analysis",
                      paste0("DEGenes_Plot_", project, "_", process, "_", automataID, ".pdf")),
           p,
           width = 6.5,
           height = 4,
           units = "in")

    DEGenes.plot[[i]] <- p
    i <- i + 1


  }

  write_csv(summary, here::here("TCGA_results", project, automataID,
                                   paste0("Summary_", project, "_", automataID, ".csv")))

  save(DEGenes.plot, file = paste0("TCGA_results/", project, "/", automataID,
                                   "/report_plots/deGenes_Plot_", project, "_", automataID, ".Rdata"))

}
