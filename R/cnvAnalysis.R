#' Generates CNV aberrations per gene
#'
#' Creates genomic ranges for the CNV data downloaded from the TCGA and finds overlaps with genomic ranges
#' originated from ENSEMBL.
#' @param project Name of the project to analyse according to the TCGA (e.g. TCGA-PRAD)
#' @param automataID Custom ID tag
#' @import magrittr
#' @import here
#' @importFrom readr read_csv write_csv
#' @importFrom dplyr filter mutate rename
#' @importFrom tidyr drop_na
#' @importFRom GenomicRanges makeGRangesFromDataFrame findOverlaps subjectHits queryHit

generateAberrations <- function(project, automataID){

  # Imports and cleans CNV data
  cnvData <- read_csv(here::here("TCGA_results", project, automataID, "TCGA_data",
                                 paste0("CNV_data_", project, "_", automataID, ".csv"))) %>%
    dplyr::select(-c(X1)) %>%
    mutate(Chromosome = gsub("X", "23", Chromosome),
           Chromosome = gsub("Y", "24", Chromosome))

  # Filters and assigns a value according to aberration type
  cnvData_filtered <- cnvData %>%
    filter(Num_Probes >= 10) %>%
    mutate(Aberration = case_when(
      Segment_Mean > 2 ~ 1,
      Segment_Mean < -2 ~ 0,
      TRUE ~ 2
    )) %>%
    filter(Aberration != 2) %>%
    dplyr::select(Sample, Chromosome, Aberration, Start, End) %>%
    drop_na()

  # Creates genomic ranges of the CNV data
  cnvGenomicRange <- makeGRangesFromDataFrame(cnvData_filtered, keep.extra.columns = TRUE)


  # Creates genomic ranges of the CNV data
  data(genesRanges)

  genesRanges.transformed <- genesRanges[genesRanges[,1]!="" & genesRanges[,2]%in%c(1:22,"X","Y"),]
  xidx <- which(genesRanges.transformed[,2]=="X")
  yidx <- which(genesRanges.transformed[,2]=="Y")
  genesRanges.transformed[xidx, 2] <- 23
  genesRanges.transformed[yidx, 2] <- 24
  genesRanges.transformed[,2] <- sapply(genesRanges.transformed[,2],as.integer)
  genesRanges.transformed <- genesRanges.transformed[order(genesRanges.transformed[,3]),]
  genesRanges.transformed <- genesRanges.transformed[order(genesRanges.transformed[,2]),]
  genesRanges.transformed  <- dplyr::rename(genesRanges.transformed, Chromosome = chromosome_name, Start = start_position,
                                            End = end_position)

  # Creates genomic ranges of the biomart data
  genesRanges.transformed <- makeGRangesFromDataFrame(genesRanges.transformed, keep.extra.columns = TRUE)

  # Find overlaps between both genomic ranges
  hits <- GenomicRanges::findOverlaps(genesRanges.transformed, cnvGenomicRange, type="within")
  cnvHits <- cbind(cnvData_filtered[subjectHits(hits),], genesRanges.transformed[queryHits(hits),])

  # Creates vectors for the aberrant region and the gene region
  AberrantRegion <- paste0(cnvHits[,2],":", cnvHits[,4],"-", cnvHits[,5])
  GeneRegion <- paste0(cnvHits[,7],":", cnvHits[,8],"-", cnvHits[,9])

  # Creates data frame with the aberrations
  geneAberrations <-  cbind(cnvHits$hgnc_symbol, cnvHits[,c(2,3)], AberrantRegion, GeneRegion, cnvHits[,c(1)]) %>%
    dplyr::rename(Sample = "cnvHits[, c(1)]", hgnc_symbol = "cnvHits$hgnc_symbol")

  geneAberrations[geneAberrations[,3]==0,3] <- "Del"
  geneAberrations[geneAberrations[,3]==1,3] <- "Amp"
  rownames(geneAberrations) <- NULL

  createTCGAfolder(project, automataID)

  write_csv(geneAberrations, here::here("TCGA_results", project, automataID, "cnv_analysis",
                                        paste0("geneAberrations_", project, "_", automataID, ".csv")))

  geneAberrations

}

# aberrations_raw <- read_csv(here::here("TCGA_results", project, automataID, "cnv_analysis",paste0("geneAberrations_", project, "_", automataID, ".csv")))

#' Performs CNV analysis
#'
#' Compares the proportion and location of chromosomical aberrations beetween each LPD groups.
#' @param project Name of the project to analyse according to the TCGA (e.g. TCGA-PRAD)
#' @param automataID Custom ID tag
#' @import magrittr
#' @importFrom dplyr mutate filter slice tally
#' @importFrom tidyr replace_na
#' @importFrom readr read_csv
#' @importFrom GenomicRanges seqinfo
#' @import BSgenome.Hsapiens.UCSC.hg38
#' @export

cnvAnalysis <- function(project, automataID){

  message("\n\n")
  message("-------------------------------------------------------")
  message("-----------------------CNV ANALYSIS--------------------")
  message("-------------------------------------------------------")
  message("\n\n")

  # Generates all the needed forlders for automata
  createTCGAfolder(project, automataID)

  # Creates a variable to store the summary or uses one already generated
  if(file.exists(here::here("TCGA_results", project, automataID,
                            paste0("Summary_", project, "_", automataID, ".csv")))){

    summary <- read_csv(here::here("TCGA_results", project, automataID,
                                   paste0("Summary_", project, "_", automataID, ".csv"))) %>%
      mutate(across(everything(), as.character))
  } else{
    summary <- tibble()
  }

  # Creates the file with aberrations per gene and sample
  aberrations_raw <- generateAberrations(project, automataID)

  # Cleans possible duplicated samples
  sampleNames <- data.frame("Barcode" = aberrations_raw$Sample) %>%
    unique() %>%
    separate(Barcode, c("Project", "TSS", "Participant", "Sample_Vial", "Portion_Analyte", "Plate", "Center"), "-", remove = FALSE) %>%
    separate(Sample_Vial, c("Sample", "Vial"), "(?<=[0-9])(?=[A-Z])") %>%
    separate(Portion_Analyte, c("Portion", "Analyte"), "(?<=[0-9])(?=[A-Z])") %>%
    group_by(TSS, Participant, Sample) %>%
    arrange(Vial, desc(Plate)) %>%
    dplyr::slice(1) %>%
    ungroup() %>%
    dplyr::select(-c(Project, TSS, Participant, Sample, Vial, Portion, Analyte, Plate, Center)) %>%
    pull()

  aberrations <- aberrations_raw %>%
    filter(Sample %in% sampleNames) %>%
    filter(!is.na(Chromosome), !is.na(hgnc_symbol)) %>%
    mutate(Sample = substr(Sample, 1, 16)) %>%
    mutate(StartPosition = sub(".*:([^.]+)-.*", "\\1", AberrantRegion)) %>%
    mutate(EndPosition = sub(".*-([^.]+).*", "\\1", AberrantRegion)) %>%
    mutate(RegionSize = as.integer(EndPosition) - as.integer(StartPosition))


  # Imports the sample list
  gammaValues <- read_csv(here::here("TCGA_results", project, automataID, "lpd_parameters",
                                     paste0("Sample_assigned_", project, "_", automataID, ".csv"))) %>%
    mutate(Sample = substr(Aliquot, 1, 16))


  # Extracts all the processes found for the dataset
  processes <- unique(gammaValues$Max_LPD)

  # Get the chromsome length to do proportions bp aberrated / bp not aberrated
  chromosomeLength <- data.frame(Length = seqlengths(Hsapiens)[1:24]) %>%
    mutate(Chromosome = rownames(.),
           Chromosome = gsub("X", "23", Chromosome),
           Chromosome = gsub("Y", "24", Chromosome),
           Chromosome = substr(Chromosome, 4, 5))

  for(process in processes){

    message("o Starting with process ", process)


    # Calculates the proportion of aberrated bp per sample and chromosome
    aberrationsProportion <- aberrations %>%
      distinct(AberrantRegion, .keep_all = TRUE) %>%
      mutate(Chromosome = as.character(Chromosome)) %>%
      group_by(Sample, Chromosome, Aberration) %>%
      summarise(Affected_bp = sum(RegionSize)) %>%
      left_join(chromosomeLength, by= "Chromosome") %>%
      mutate(ProportionAffected = Affected_bp/Length) %>%
      mutate(Percentage = ProportionAffected * 100) %>%
      mutate(SelectedGroup = case_when(
        Sample %in% filter(gammaValues, Max_LPD == process)$Sample ~ 1,
        TRUE ~ 0)) %>%
      arrange(SelectedGroup, Aberration, Sample, Chromosome)

    if(length(unique(filter(aberrationsProportion, Aberration == "Amp")$SelectedGroup))  == 2){
      aberrationsAmp.p <- wilcox.test(formula = ProportionAffected ~ as.factor(SelectedGroup),
                                      data = filter(aberrationsProportion, Aberration == "Amp"))$p.value
    } else{
      aberrationsAmp.p <- NA
    }

    if(length(unique(filter(aberrationsProportion, Aberration == "Del")$SelectedGroup))  == 2){
    aberrationsDel.p <- wilcox.test(formula = ProportionAffected ~ as.factor(SelectedGroup),
                                    data = filter(aberrationsProportion, Aberration == "Del"))$p.value
    } else{
      aberrationsDel.p <- NA
    }

    write_csv(aberrationsProportion, here::here("TCGA_results", project, automataID, "cnv_analysis",
                                                paste0("aberrationsProportion_", project, "_",
                                                       process, "_", automataID, ".csv")))



    # Calculates proportion of aberration for each sample ----------------------
    aberrationsSample <- aberrationsProportion %>%
      group_by(Sample, Aberration, SelectedGroup) %>%
      summarise(Affected_bp.sum = sum(Affected_bp, na.rm = TRUE), Length.sum = sum(as.numeric(Length), na.rm = TRUE),
                ProportionAffected.total = Affected_bp.sum/Length.sum) %>%
      arrange(SelectedGroup, Aberration, Sample)

    if(length(unique(filter(aberrationsSample, Aberration == "Amp")$SelectedGroup)) == 2){
      aberrationsSampleAmp.p <- wilcox.test(formula = ProportionAffected.total ~ as.factor(SelectedGroup),
                                            data = filter(aberrationsSample, Aberration == "Amp"))$p.value
    } else {
      aberrationsSampleAmp.p <- NA
    }

    if(length(unique(filter(aberrationsSample, Aberration == "Del")$SelectedGroup)) == 2){
    aberrationsSampleDel.p <- wilcox.test(formula = ProportionAffected.total ~ as.factor(SelectedGroup),
                                          data = filter(aberrationsSample, Aberration == "Del"))$p.value
    } else {
      aberrationsSampleDel.p <- NA
    }

    write_csv(aberrationsSample, here::here("TCGA_results", project, automataID, "cnv_analysis",
                                            paste0("aberrationsSample_", project, "_", process, "_", automataID, ".csv")))



    # Plots it
    aberrationSample.plotData <- dplyr::mutate(aberrationsSample, SelectedGroup = as.factor(ifelse(SelectedGroup == 1, process, "Other_LPD")))

    plot <- ggplot(aberrationSample.plotData,
                   aes(x = paste(Aberration, SelectedGroup, sep = "+"), y = ProportionAffected.total, fill = SelectedGroup)) +
      geom_boxplot() +
      ylab("Proportion of bp affected by aberrations") +
      coord_flip() +
      theme_minimal() +
      scale_fill_manual(values = c("darkslategray4", 'aquamarine3')) +
      ggtitle(paste("Affected bp for", process, project),
              subtitle = paste0("p-value Amp: ", round(aberrationsSampleAmp.p, digits = 5), "\n",
                                "p-value Del: ", round(aberrationsSampleDel.p, digits = 5))) +
      ylim(0, 1) +
      theme(axis.title.y = element_blank(),
            plot.title = element_text(hjust = 0.5, face = "bold"),
            legend.position = "none")

    ggsave(here::here("TCGA_results", project, automataID, "cnv_analysis",
                      paste0("aberrationSample_Plot_", project, "_", process, "_", automataID, ".pdf")),
           plot, width = 6, height = 4, units = "in")

    aberrationSample_plot <- plot

    save(aberrationSample_plot, file = paste0("TCGA_results/", project, "/", automataID,
                                              "/report_plots/aberrationSample_plot_", project, "_", process, "_", automataID, ".Rdata"))


    # Calculates proportion of aberration for each chromosome ----------------------

    # Calculates the proportion per chromosome for each sample
    aberrationChromosome <- aberrationsProportion %>%
      group_by(Sample, Chromosome, Aberration, SelectedGroup) %>%
      summarise(Affected_bp.sum = sum(as.numeric(Affected_bp), na.rm = TRUE), Length.sum = sum(as.numeric(Length), na.rm = TRUE),
                ProportionAffected.total = Affected_bp.sum/Length.sum) %>%
      dplyr::mutate(SelectedGroup = as.factor(ifelse(SelectedGroup == 1, process, "Other_LPD"))) %>%
      ungroup()

    write_csv(aberrationChromosome, here::here("TCGA_results", project, automataID, "cnv_analysis",
                                               paste0("aberrationsChromosome_", project, "_", process, "_", automataID, ".csv")))

    # For amplification only ######
    aberrationChromosome.amp <- filter(aberrationChromosome, Aberration == "Amp")

    # Calculates if there are significant differences
    chromo.pvalues <- tibble()

    for(i in as.numeric(unique(aberrationChromosome.amp$Chromosome))){
      chromo.stats <- filter(aberrationChromosome.amp, Chromosome == i)
      if(length(unique(chromo.stats$SelectedGroup)) != 2){
        next()
      }
      chromo.wilcox <- wilcox.test(formula = ProportionAffected.total ~ as.factor(SelectedGroup), data=chromo.stats)$p.value
      chromo.pvalues <- bind_rows(chromo.pvalues, c("Chromosome" = as.character(i), "P.val" = chromo.wilcox))
    }

    # Assigns asterisks when they is lower than 0.05 and duplicates to make it compatible with ggplot

    if(nrow(chromo.pvalues) > 0){

      chromo.pvalues <- chromo.pvalues %>%
        arrange(as.integer(Chromosome)) %>%
        mutate(Significant = ifelse(P.val <= 0.05, "*", ""),
               SelectedGroup = process)

      chromoClon <- chromo.pvalues %>%
        mutate(Significant = "",
               SelectedGroup = "Other_LPD")

      chromo.pvalues <- bind_rows(chromo.pvalues, chromoClon) %>%
        mutate(y.position = 1)

      # Plots it
      plot <- ggplot(aberrationChromosome.amp,
                     aes(x  = reorder(Chromosome, as.integer(Chromosome)), y = ProportionAffected.total, fill = SelectedGroup)) +
        geom_boxplot() +
        geom_text(data = chromo.pvalues, aes(label = Significant, y = y.position, colour = "red"), size = 12, show.legend = FALSE) +
        theme_minimal() +
        ggtitle(paste("Amplified bp for ", process, project)) +
        ylab("Proportion of bps affected by amplification") +
        xlab("Chromosome") +
        ylim(0, 1) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
        scale_fill_manual(values=c("darkslategray4", 'aquamarine3'), labels = c(process, "Other LPD"))

      ggsave(here::here("TCGA_results", project, automataID, "cnv_analysis",
                        paste0("aberrationChromosomePlot_AMP_", project, "_", process, "_", automataID, ".pdf")),
             plot, width = 12, height = 6, units = "in")

      aberrationChromosomePlot_Amp <- plot

      save(aberrationChromosomePlot_Amp, file = paste0("TCGA_results/", project, "/", automataID,
                                                       "/report_plots/aberrationChromosome_Amp_plot_", project, "_", process, "_", automataID, ".Rdata"))


    }


    # For deletion only ####
    aberrationChromosome.del <- filter(aberrationChromosome, Aberration == "Del")

    # Calculates if there are significant differences
    chromo.pvalues <- tibble()

    for(i in as.numeric(unique(aberrationChromosome.del$Chromosome))){
      chromo.stats <- filter(aberrationChromosome.del, Chromosome == i)
      if(length(unique(chromo.stats$SelectedGroup)) < 2){
        next()
      }
      chromo.wilcox <- wilcox.test(formula = ProportionAffected.total ~ as.factor(SelectedGroup), data=chromo.stats)$p.value
      chromo.pvalues <- bind_rows(chromo.pvalues, c("Chromosome" = as.character(i), "P.val" = chromo.wilcox))
    }

    # Assigns asterisks when they is lower than 0.05 and duplicates to make it compatible with ggplot
    if(nrow(chromo.pvalues) > 0){

      chromo.pvalues <- chromo.pvalues %>%
        arrange(as.integer(Chromosome)) %>%
        mutate(Significant = ifelse(P.val <= 0.05, "*", ""),
               SelectedGroup = process)

      chromoClon <- chromo.pvalues %>%
        mutate(Significant = "",
               SelectedGroup = "Other_LPD")

      chromo.pvalues <- bind_rows(chromo.pvalues, chromoClon) %>%
        mutate(y.position = 1)

      # Plots it
      plot <- ggplot(aberrationChromosome.del,
                     aes(x  = reorder(Chromosome, as.integer(Chromosome)), y = ProportionAffected.total, fill = SelectedGroup)) +
        geom_boxplot() +
        geom_text(data = chromo.pvalues, aes(label = Significant, y = y.position, colour = "red"), size = 12, show.legend = FALSE) +
        theme_minimal() +
        ggtitle(paste("Deleted bp for ", process, project)) +
        ylab("Proportion of bps affected by deletion") +
        xlab("Chromosome") +
        ylim(0, 1) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
        scale_fill_manual(values=c("darkslategray4", 'aquamarine3'), labels = c(process, "Other LPD"))

      ggsave(here::here("TCGA_results", project, automataID, "cnv_analysis",
                        paste0("aberrationChromosomePlot_DEL_", project, "_", process, "_", automataID, ".pdf")),
             plot, width = 12, height = 6, units = "in")

      aberrationChromosomePlot_Del <- plot

      save(aberrationChromosomePlot_Del, file = paste0("TCGA_results/", project, "/", automataID,
                                                       "/report_plots/aberrationChromosome_Del_plot_", project, "_", process, "_", automataID, ".Rdata"))

    }



    # Extracts which genes are affected and runs pathway analysis
    aberratedGenes <- aberrations %>%
      group_by(hgnc_symbol, Aberration, Sample) %>%
      dplyr::slice(1) %>%
      ungroup() %>%
      right_join(gammaValues)

    aberratedGenes.Selected <- aberratedGenes %>%
      filter(Max_LPD == process) %>%
      dplyr::select(hgnc_symbol, Aberration, Sample)

    aberratedGenes.Context <- aberratedGenes %>%
      filter(Max_LPD != process) %>%
      dplyr::select(hgnc_symbol, Aberration, Sample)

    aberratedGenes_amp <- filter(aberratedGenes.Selected, Aberration == "Amp")
    if(nrow(aberratedGenes_amp) > 0){

      runPathwayAnalysis(project, automataID, process, geneNames = aberratedGenes_amp$hgnc_symbol,
                         contextNames = aberratedGenes.Selected$hgnc_symbol, positive.l2FC = TRUE, type = "cnv")
    }

    aberratedGenes_del <- filter(aberratedGenes.Selected, Aberration == "Del")
    if(nrow(aberratedGenes_del) > 0){

      runPathwayAnalysis(project, automataID, process, geneNames = aberratedGenes_del$hgnc_symbol,
                         contextNames = aberratedGenes.Selected$hgnc_symbol, positive.l2FC = FALSE, type = "cnv")
    }

    aberratedGenesCount.Selected <- aberratedGenes.Selected %>%
      group_by(hgnc_symbol, Aberration) %>%
      tally() %>%
      spread(key = Aberration, value = n) %>%
      tidyr::replace_na(list(hgnc_symbol = "Unknown", Amp = 0, Del = 0))

    if(NA  %in% colnames(aberratedGenesCount.Selected)){
      aberratedGenesCount.Selected <- dplyr::select(aberratedGenesCount.Selected, -`<NA>`)
    }

    if(ncol(aberratedGenesCount.Selected) != 3){
      next()
    }

    aberratedGenesCount.Context <- aberratedGenes.Context %>%
      group_by(hgnc_symbol, Aberration) %>%
      tally() %>%
      spread(key = Aberration, value = n) %>%
      tidyr::replace_na(list(hgnc_symbol = "Unknown", Amp = 0, Del = 0))

    if(NA  %in% colnames(aberratedGenesCount.Context)){
      aberratedGenesCount.Context<- dplyr::select(aberratedGenesCount.Context, -`<NA>`)
    }

    if(ncol(aberratedGenesCount.Context) != 3){
      next()
    }

    numberSamples.Selected <- length(unique(aberratedGenes.Selected$Sample))
    numberSamples.Context <- length(unique(aberratedGenes.Context$Sample))

    aberratedComparison <- aberratedGenesCount.Selected %>%
      left_join(aberratedGenesCount.Context, by = "hgnc_symbol") %>%
      dplyr::rename(Amp.Selected = "Amp.x", Del.Selected = "Del.x", Amp.Context = "Amp.y", Del.Context = "Del.y") %>%
      mutate(notAmp.Selected = numberSamples.Selected - Amp.Selected,
             notAmp.Context = numberSamples.Context - Amp.Context,
             notDel.Selected = numberSamples.Selected - Del.Selected,
             notDel.Context = numberSamples.Context - Del.Context) %>%
      dplyr::select(hgnc_symbol, Amp.Selected, Amp.Context, notAmp.Selected,  notAmp.Context, everything()) %>%
      dplyr::filter(!is.na(hgnc_symbol)) %>%
      replace(is.na(.), 0)

    aberratedComparison$Amp.chisq <- apply(aberratedComparison, 1, function(x){
      tbl = matrix(as.numeric(x[2:5]), ncol = 2, byrow = TRUE)
      chisq.test(tbl)$p.value
    })

    aberratedComparison$Amp.padj <- p.adjust(aberratedComparison$Amp.chisq, method = "BH")

    aberratedComparison$Del.chisq <- apply(aberratedComparison, 1, function(x){
      tbl = matrix(as.numeric(x[6:9]), ncol = 2, byrow = TRUE)
      chisq.test(tbl)$p.value
    })

    aberratedComparison$Del.padj <- p.adjust(aberratedComparison$Del.chisq, method = "BH")

    write_csv(aberratedComparison, here::here("TCGA_results", project, automataID, "cnv_analysis",
                                         paste0("aberrationsGene_", project, "_", process, "_", automataID, ".csv")))

    # Creates the summary vector
    significantCount <- aberratedComparison %>%
      dplyr::filter(Amp.padj <= 0.05 | Del.padj <= 0.05)

    summaryVector <- c(Analysis = "Copy Number Variations",
                       Test = "Proportion of genes affected",
                       LPD_Group = process,
                       Significant = ifelse(nrow(significantCount) > 0, TRUE, FALSE),
                       n = nrow(significantCount),
                       p.value = NA,
                       residual = NA,
                       date = as.character(Sys.Date()))

    summary <- bind_rows(summary, summaryVector)

  }


  message("\n\n")
  message("=======================================================")
  message("\n\n")


}
