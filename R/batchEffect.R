#' Analyse possible batch effect
#'
#' Compares proportion of different TSS and sample type for each LPD group to see if there
#' is a significant diffence between them.
#'
#' @param project Name of the project to analyse according to the TCGA (e.g. TCGA-PRAD)
#' @param automataID Custom ID tag
#' @return A CSV file with proportions of TSS; a CSV filr with proportions of sample type: a
#' CSV file with proportions of TSS and sample type
#' @importFrom readr read_csv
#' @importFrom dplyr mutate group_by summarise bind_rows ungroup
#' @importFrom tidyr spread
#' @importFrom tibble remove_rownames column_to_rownames
#' @import magrittr
#' @import here
#' @import chisq.posthoc.test
#' @export


batchEffect <- function(project, automataID){

  message("\n\n")
  message("-------------------------------------------------------")
  message("-----------------------BATCH EFFECT--------------------")
  message("-------------------------------------------------------")
  message("\n\n")

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

  # Read the gamma values and creates columns with the TSS and the sample type
  gammaValues <- read_csv(here::here("TCGA_results", project, automataID, "lpd_parameters",
                                     paste0("Sample_assigned_", project, "_", automataID, ".csv"))) %>%
    mutate(TSS = substr(Sample, 6, 7),
           SampleType = substr(Sample, 14, 16))


  # Check by TSS -----------------------------
  message("o Batch effect by TSS")
  assignedTSS.proportion <- gammaValues %>%
    group_by(Max_LPD, TSS) %>%
    summarise(n = n()) %>%
    spread(key = TSS, value = n) %>%
    replace(is.na(.), 0)

  # I need to calculate the chisq of one group against all other together
  chi.pvals <- vector()
  chi.residual <- vector()

  for(i in 1:nrow(assignedTSS.proportion)){
    message("oo Process LPD_", i)
    tempMatrix <- as.matrix(bind_rows(assignedTSS.proportion[i, -1], colSums(assignedTSS.proportion[-i, -1])))

    if(ncol(tempMatrix) < 2){
      next
    }

    chisq <- chisq.test(tempMatrix)

    dimnames(tempMatrix) <- list(Process = c("LPD_Selected", "LPD_Context"),
                                 TSS = colnames(tempMatrix))

    # Uses posthoc to calculate the p.values of each column
    chis.posthoc <- chisq.posthoc.test(tempMatrix) %>%
      slice(1:2) %>%
      select(-Dimension) %>%
      remove_rownames() %>%
      column_to_rownames("Value") %>%
      t() %>%
      as.data.frame() %>%
      rownames_to_column("TSS") %>%
      filter(`p values` <= 0.05)

    # Gets only the names of the significant TSS
    significantResiduals <- ifelse(nrow(chis.posthoc) > 0, paste(chis.posthoc$TSS, collapse = ", "),
                                   "None")

    chi.pvals[i] <- chisq$p.value
    # Takes the name of the column with the highest sum of residuals as the most variant
    chi.residual[i] <- significantResiduals

  }

  if(length(chi.pvals) > 0){

    assignedTSS <- assignedTSS.proportion %>%
      ungroup() %>%
      mutate(chisq.pval = chi.pvals,
             chisq.padj = p.adjust(chisq.pval, method = "BH"),
             significant.residual = chi.residual)

    write_csv(assignedTSS, here::here("TCGA_results", project, automataID, "batch_effect",
                                      paste0("TSSProportions_", project, "_", automataID, ".csv")))


    # Creates the summary attachment
    summaryTibble <- assignedTSS %>%
      dplyr::select(Max_LPD, chisq.padj, significant.residual) %>%
      rename(LPD_Group = "Max_LPD", p.value = "chisq.padj", residual = "significant.residual") %>%
      mutate(Analysis = "Batch effect",
             Test = "TSS proportion",
             Significant = ifelse(p.value <= 0.05, TRUE, FALSE),
             n = NA,
             date = as.character(Sys.Date())) %>%
      dplyr::select(Analysis, Test, LPD_Group, Significant, n, p.value, residual, date) %>%
      mutate(across(everything(), as.character))

    summary <- bind_rows(summary, summaryTibble)

  }

  # Check by sample type -----------------------------
  assignedsampleType.proportion <- gammaValues %>%
    group_by(Max_LPD, SampleType) %>%
    summarise(n = n()) %>%
    spread(key = SampleType, value = n) %>%
    replace(is.na(.), 0)

  # I need to calculate the chisq of one group against all other together
  chi.pvals <- vector()
  chi.residual <- vector()

  for(i in 1:nrow(assignedsampleType.proportion)){
    message("oo Process LPD_", i)
    tempMatrix <- as.matrix(bind_rows(assignedsampleType.proportion[i, -1], colSums(assignedsampleType.proportion[-i, -1])))
    chisq <- chisq.test(tempMatrix)

    dimnames(tempMatrix) <- list(Process = c("LPD_Selected", "LPD_Context"),
                                 TSS = colnames(tempMatrix))

    if(ncol(tempMatrix) == 1){
      next()
    }

    # Uses posthoc to calculate the p.values of each column
    chis.posthoc <- chisq.posthoc.test(tempMatrix) %>%
      slice(1:2) %>%
      select(-Dimension) %>%
      remove_rownames() %>%
      column_to_rownames("Value") %>%
      t() %>%
      as.data.frame() %>%
      rownames_to_column("TSS") %>%
      filter(`p values` <= 0.05)

    # Gets only the names of the significant TSS
    significantResiduals <- ifelse(nrow(chis.posthoc) > 0, paste(chis.posthoc$TSS, collapse = ", "),
                                   "None")

    chi.pvals[i] <- chisq$p.value
    # Takes the name of the column with the highest sum of residuals as the most variant
    chi.residual[i] <- significantResiduals
  }

  if(ncol(assignedsampleType.proportion) > 2){

    assignedsampleType <- assignedsampleType.proportion %>%
      ungroup() %>%
      mutate(chisq.pval = chi.pvals,
             chisq.padj = p.adjust(chisq.pval, method = "BH"),
             significant.residual = chi.residual)

    write_csv(assignedsampleType, here::here("TCGA_results", project, automataID, "batch_effect",
                                             paste0("sampleTypeProportions_", project, "_", automataID, ".csv")))


    # Creates the summary attachment
    summaryTibble <- assignedsampleType %>%
      dplyr::select(Max_LPD, chisq.padj, significant.residual) %>%
      rename(LPD_Group = "Max_LPD", p.value = "chisq.padj", residual = "significant.residual") %>%
      mutate(Analysis = "Batch effect",
             Test = "Sample type proportion",
             Significant = ifelse(p.value <= 0.05, TRUE, FALSE),
             n = NA,
             date = as.character(Sys.Date())) %>%
      dplyr::select(Analysis, Test, LPD_Group, Significant, n, p.value, residual, date) %>%
      mutate(across(everything(), as.character))

    summary <- bind_rows(summary, summaryTibble)

    }


  write_csv(summary, here::here("TCGA_results", project, automataID,
                                paste0("Summary_", project, "_", automataID, ".csv")))

message("\n\n")
message("=======================================================")
message("\n\n")



}


