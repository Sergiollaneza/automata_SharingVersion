#' Clinical analysis
#'
#' Does a clinical analysis
#'
#' @param .
#' @param .
#' @importFrom readr read_csv
#' @importFrom dplyr filter mutate inner_join
#' @importFrom tidyr spread
#' @importFrom survminer ggsurvplot surv_pvalue
#' @importFrom survival survfit Surv
#' @importFrom tibble remove_rownames rownames_to_column
#' @import here
#' @import chisq.posthoc.test
#' @export

clinicalAnalysis <- function(project, automataID){

  message("\n\n")
  message("-------------------------------------------------------")
  message("------------------CLINICAL ANALYSIS--------------------")
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

  # Read the data
  clinical <- read_csv(here::here("TCGA_results", project, automataID, "TCGA_data",
                                  paste0("clinical_data_", project, "_", automataID, ".csv")))

  gammaValues <- read_csv(here::here("TCGA_results", project, automataID, "lpd_parameters",
                                     paste0("Sample_assigned_", project, "_", automataID, ".csv"))) %>%
    filter(substr(Sample, 14, 15) != "11") %>%
    mutate(Sample = substr(Sample, 1,12))

  clinicalAssigned <- inner_join(clinical, gammaValues, by = c("bcr_patient_barcode" = "Sample"))


  # Proportion analysis -----------------------------------------------------

  message("o Proportion analysis")

  # By cancer stage -----------------------------

  message("oo By cancer stage")

  if(!(project %in% c("TCGA-LGG", "TCGA-GBM", "TCGA-LAML", "TCGA-SARC", "TCGA-PCPG"))){

    pathologicStage.proportion <- clinicalAssigned %>%
      group_by(Max_LPD, stage_event_pathologic_stage) %>%
      summarise(n = n()) %>%
      filter(is.na(stage_event_pathologic_stage) == FALSE) %>%
      spread(key = stage_event_pathologic_stage, value = n) %>%
      replace(is.na(.), 0)

    if(nrow(pathologicStage.proportion) > 0){

      # I need to calculate the chisq of one group against all other together
      chi.pvals <- vector()
      chi.residual <- vector()

      for(i in 1:nrow(pathologicStage.proportion)){

        tempMatrix <- as.matrix(bind_rows(pathologicStage.proportion[i, -1], colSums(pathologicStage.proportion[-i, -1])))
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

      pathologicalStage <- pathologicStage.proportion %>%
        ungroup() %>%
        mutate(chisq.pval = chi.pvals,
               chisq.padj = p.adjust(chisq.pval, method = "BH"),
               highest.residual = chi.residual)

      write_csv(pathologicalStage, here::here("TCGA_results", project, automataID, "clinical_analysis",
                                              paste0("pathologicalStageProportions_", project, "_", automataID, ".csv")))

      # Creates the summary attachment
      summaryTibble <- pathologicalStage %>%
        dplyr::select(Max_LPD, chisq.padj) %>%
        rename(LPD_Group = "Max_LPD", p.value = "chisq.padj") %>%
        mutate(Analysis = "Clinical analysis",
               Test = "Pathological Stage",
               Significant = ifelse(p.value <= 0.05, TRUE, FALSE),
               n = NA,
               date = as.character(Sys.Date()),
               residual = NA) %>%
        dplyr::select(Analysis, Test, LPD_Group, Significant, n, p.value, residual, date) %>%
        mutate(across(everything(), as.character))

      summary <- bind_rows(summary, summaryTibble)

    }
  } else{
    message("ooo Not available for this cancer")
  }



  # By gender -----------------------------

  message("oo By gender")

  gender.proportion <- clinicalAssigned %>%
    group_by(Max_LPD, gender) %>%
    summarise(n = n()) %>%
    filter(is.na(gender) == FALSE) %>%
    spread(key = gender, value = n) %>%
    replace(is.na(.), 0)

  if(ncol(gender.proportion) > 2){

    # I need to calculate the chisq of one group against all other together
    chi.pvals <- vector()
    chi.residual <- vector()

    for(i in 1:nrow(gender.proportion)){

      tempMatrix <- as.matrix(bind_rows(gender.proportion[i, -1], colSums(gender.proportion[-i, -1])))
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

    gender <- gender.proportion %>%
      ungroup() %>%
      mutate(chisq.pval = chi.pvals,
             chisq.padj = p.adjust(chisq.pval, method = "BH"),
             highest.residual = chi.residual)

    write_csv(gender, here::here("TCGA_results", project, automataID, "clinical_analysis",
                                 paste0("genderProportions_", project, "_", automataID, ".csv")))


    # Creates the summary attachment
    summaryTibble <- gender %>%
      dplyr::select(Max_LPD, chisq.padj, highest.residual) %>%
      rename(LPD_Group = "Max_LPD", p.value = "chisq.padj", residual = "highest.residual") %>%
      mutate(Analysis = "Clinical Analysis",
             Test = "Gender proportion",
             Significant = ifelse(p.value <= 0.05, TRUE, FALSE),
             n = NA,
             date = as.character(Sys.Date())) %>%
      dplyr::select(Analysis, Test, LPD_Group, Significant, n, p.value, residual, date) %>%
      mutate(across(everything(), as.character))

    summary <- bind_rows(summary, summaryTibble)


  }

  # By race -----------------------------

  message("oo By race")

  race.proportion <- clinicalAssigned %>%
    group_by(Max_LPD, race_list) %>%
    summarise(n = n()) %>%
    filter(is.na(race_list) == FALSE) %>%
    spread(key = race_list, value = n) %>%
    replace(is.na(.), 0)

  # I need to calculate the chisq of one group against all other together
  chi.pvals <- vector()
  chi.residual <- vector()

  for(i in 1:nrow(race.proportion)){

    tempMatrix <- as.matrix(bind_rows(race.proportion[i, -1], colSums(race.proportion[-i, -1])))
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

  race <- race.proportion %>%
    ungroup() %>%
    mutate(chisq.pval = chi.pvals,
           chisq.padj = p.adjust(chisq.pval, method = "BH"),
           highest.residual = chi.residual)

  write_csv(race, here::here("TCGA_results", project, automataID, "clinical_analysis",
                             paste0("raceProportions_", project, "_", automataID, ".csv")))

  # Creates the summary attachment
  summaryTibble <- race %>%
    dplyr::select(Max_LPD, chisq.padj, highest.residual) %>%
    rename(LPD_Group = "Max_LPD", p.value = "chisq.padj", residual = "highest.residual") %>%
    mutate(Analysis = "Clinical Analysis",
           Test = "Race proportion",
           Significant = ifelse(p.value <= 0.05, TRUE, FALSE),
           n = NA,
           date = as.character(Sys.Date())) %>%
    dplyr::select(Analysis, Test, LPD_Group, Significant, n, p.value, residual, date) %>%
    mutate(across(everything(), as.character))

  summary <- bind_rows(summary, summaryTibble)


  # By primary pathology -----------------------------
  if(project == "TCGA-CHOL"){

    message("oo By primary pathology")

    pathology.proportion <- clinicalAssigned %>%
      group_by(Max_LPD, primary_pathology_histological_type) %>%
      summarise(n = n()) %>%
      filter(is.na(primary_pathology_histological_type) == FALSE) %>%
      spread(key = primary_pathology_histological_type, value = n) %>%
      replace(is.na(.), 0)

    # I need to calculate the chisq of one group against all other together
    chi.pvals <- vector()
    chi.residual <- vector()

    for(i in 1:nrow(pathology.proportion)){

      tempMatrix <- as.matrix(bind_rows(pathology.proportion[i, -1], colSums(pathology.proportion[-i, -1])))
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

    pathology <- pathology.proportion %>%
      ungroup() %>%
      mutate(chisq.pval = chi.pvals,
             chisq.padj = p.adjust(chisq.pval, method = "BH"),
             highest.residual = chi.residual)

    write_csv(pathology, here::here("TCGA_results", project, automataID, "clinical_analysis",
                                    paste0("pathologyProportions_", project, "_", automataID, ".csv")))

  }

  # By gleason score for TCGA-PRAD  -----------------------------

  message("oo By gleason score")

  if(project == "TCGA-PRAD"){

    gleason.proportion <- clinicalAssigned %>%
      mutate(stage_event_gleason_grading = substr(stage_event_gleason_grading, 0, 1),
             stage_event_gleason_grading = ifelse(stage_event_gleason_grading == 1,
                                                  10, stage_event_gleason_grading)) %>%
      group_by(Max_LPD, stage_event_gleason_grading) %>%
      summarise(n = n()) %>%
      filter(is.na(stage_event_gleason_grading) == FALSE) %>%
      spread(key = stage_event_gleason_grading, value = n) %>%
      replace(is.na(.), 0)

    # I need to calculate the chisq of one group against all other together
    chi.pvals <- vector()
    chi.residual <- vector()

    for(i in 1:nrow(gleason.proportion)){

      tempMatrix <- as.matrix(bind_rows(gleason.proportion[i, -1], colSums(gleason.proportion[-i, -1])))
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

    gleason <- gleason.proportion %>%
      ungroup() %>%
      mutate(chisq.pval = chi.pvals,
             chisq.padj = p.adjust(chisq.pval, method = "BH"),
             highest.residual = chi.residual)

    write_csv(gleason, here::here("TCGA_results", project, automataID, "clinical_analysis",
                                  paste0("gleasonProportions_", project, "_", automataID, ".csv")))

  } else {
    message("ooo Not available for this cancer")
  }


  # Survival analysis -----------------------------------------------------

  message("o Survival analysis")

  # Does the survival analysis overall
  if(project == "TCGA-PRAD"){

    message("ooo Project is prostate cancer")

    survivalData <- data.frame(Sample_ID = clinicalAssigned$bcr_patient_barcode,
                               LPD_Group = clinicalAssigned$Max_LPD,
                               Event = clinicalAssigned$biochemical_recurrence,
                               Time = clinicalAssigned$days_to_last_followup,
                               BCR = clinicalAssigned$days_to_first_biochemical_recurrence) %>%
      filter(is.na(Event) == FALSE, is.na(Time) == FALSE)

    # Change values to BCR
    survivalData$Time <- ifelse(survivalData$Event == "NO", survivalData$Time, survivalData$BCR)
    survivalData$BCR  = NULL

    # Change Yes and No to 1 and 0
    survivalData$Event = ifelse(survivalData$Event == "YES", 1, 0)

    # Keeping the LPD groups order, sorts the patients by hours (from less to more)
    survivalData <- survivalData[order(survivalData$LPD_Group, survivalData$Time),]

    # Create the km estimations
    km_fit <- survminer::surv_fit(Surv(time = Time/360, event = Event) ~ LPD_Group, data = survivalData)

  } else {

    message("ooo Project is not prostate cancer")

    # The days to death will be used as event

    message("ooo Setting days to death as event")

    survivalData <- data.frame(Sample_ID = clinicalAssigned$bcr_patient_barcode,
                               LPD_Group = clinicalAssigned$Max_LPD,
                               Time = clinicalAssigned$days_to_death,
                               FollowUp = clinicalAssigned$days_to_last_followup,
                               Event = clinicalAssigned$vital_status) %>%
      mutate(Time = ifelse(is.na(Time), FollowUp, Time)) %>%
      dplyr::filter(is.na(Event) == FALSE, is.na(Time) == FALSE) %>%
      dplyr::select(-FollowUp) %>%
      mutate(Event = ifelse(Event == "Dead", TRUE, FALSE))


    # Create the km estimations

    message("ooo Generating KM estimations")

    km_fit <- survminer::surv_fit(Surv(time = Time/360, event = Event) ~ LPD_Group, data = survivalData)
  }


  survivalOverall_plot <- ggsurvplot(km_fit, pval = TRUE, xlab = "Years", title = paste0("Kaplan-Meier estimation for all signatures ", project, " ", automataID),
                                     legend = "bottom", risk.table = TRUE)

  pdf(here::here("TCGA_results", project, automataID, "clinical_analysis",
                 paste0("KM_Survival_Overall_Plot_",  automataID, ".pdf")),
      width = 10, height = 10, onefile = FALSE)

  print(survivalOverall_plot)
  dev.off()

  save(survivalOverall_plot, file = paste0("TCGA_results/", project, "/", automataID,
                                           "/report_plots/survivalOverall_plot_",
                                           project, "_", automataID, ".Rdata"))


  # Creates the summary attachment
  fit_pval <- surv_pvalue(km_fit)$pval

  summaryVector <- c(Analysis = "Clinical Analysis",
                     Test = "Survival Overall",
                     LPD_Group = "All",
                     Significant = ifelse(fit_pval <= 0.05, TRUE, FALSE),
                     n = NA,
                     p.value = fit_pval,
                     residual = NA,
                     date = as.character(Sys.Date()))


  summary <- bind_rows(summary, summaryVector)


  # Now survival plot for one vs all ----------------------

  # Extracts all the processes found for the dataset
  processes <- unique(gammaValues$Max_LPD)
  chi.pvalues <- tibble()
  wilcox.pvalues <- tibble()

  for(process in processes){

    clinicalAssigned_bin <- dplyr::mutate(clinicalAssigned,
                                          SelectedGroup = ifelse(Max_LPD == process, process, "Other_LPD"))



    # Analysis for prostate
    if(project == "TCGA-PRAD"){

      survivalData <- data.frame(Sample_ID = clinicalAssigned_bin$bcr_patient_barcode,
                                 LPD_Group = clinicalAssigned_bin$Max_LPD,
                                 Event = clinicalAssigned_bin$biochemical_recurrence,
                                 Time = clinicalAssigned_bin$days_to_last_followup,
                                 BCR = clinicalAssigned_bin$days_to_first_biochemical_recurrence,
                                 SelectedGroup = clinicalAssigned_bin$SelectedGroup) %>%
        filter(is.na(Event) == FALSE, is.na(Time) == FALSE)

      # Change values to BCR
      survivalData$Time <- ifelse(survivalData$Event == "NO", survivalData$Time, survivalData$BCR)
      survivalData$BCR  = NULL

      # Change Yes and No to 1 and 0
      survivalData$Event = ifelse(survivalData$Event == "YES", 1, 0)

      # Keeping the LPD groups order, sorts the patients by hours (from less to more)
      survivalData <- survivalData[order(survivalData$SelectedGroup, survivalData$Time),]

      # Create the km estimations
      km_fit <- survminer::surv_fit(Surv(time = Time/360, event = Event) ~ SelectedGroup, data = survivalData)

    } else {

      # The days to death will be used as event
      survivalData <- data.frame(Sample_ID = clinicalAssigned_bin$bcr_patient_barcode,
                                 LPD_Group = clinicalAssigned_bin$Max_LPD,
                                 Time = clinicalAssigned_bin$days_to_death,
                                 FollowUp = clinicalAssigned_bin$days_to_last_followup,
                                 Event = clinicalAssigned_bin$vital_status,
                                 SelectedGroup = clinicalAssigned_bin$SelectedGroup) %>%
        mutate(Time = ifelse(is.na(Time), FollowUp, Time)) %>%
        dplyr::filter(is.na(Event) == FALSE, is.na(Time) == FALSE) %>%
        dplyr::select(-FollowUp) %>%
        mutate(Event = ifelse(Event == "Dead", TRUE, FALSE))


      # Create the km estimations
      km_fit <- survminer::surv_fit(Surv(time = Time/360, event = Event) ~ SelectedGroup, data = survivalData)
    }


    survivalOverall_plot <- ggsurvplot(km_fit, pval = TRUE, xlab = "Years",
                                       title = paste0("Kaplan-Meier estimaton for ", project, " ", process, " ", automataID),
                                       legend = "bottom", risk.table = TRUE)

    pdf(here::here("TCGA_results", project, automataID, "clinical_analysis",
                   paste0("KM_Survival_Overall_Plot_",  process, "_", automataID, ".pdf")),
        width = 10, height = 10, onefile = FALSE)

    print(survivalOverall_plot)
    dev.off()

    save(survivalOverall_plot, file = paste0("TCGA_results/", project, "/", automataID,
                                             "/report_plots/survivalOverall_plot_",
                                             project, "_", process, "_", automataID, ".Rdata"))


    # Creates the summary attachment
    fit_pval <- surv_pvalue(km_fit)$pval

    summaryVector <- c(Analysis = "Clinical Analysis",
                       Test = "Survival Analyssis",
                       LPD_Group = process,
                       Significant = ifelse(fit_pval <= 0.05, TRUE, FALSE),
                       n = NA,
                       p.value = fit_pval,
                       residual = NA,
                       date = as.character(Sys.Date()))


    summary <- bind_rows(summary, summaryVector)


    # Performs extra clinical stuff for some cancers
    extraList <- extraClinical(project = project, automataID = automataID, clinicalData = clinicalAssigned_bin, process = process)

    if(!is.null(extraList)){
      # extracts the wilcox for PSA differences
      wilcoxVector <- c(Process = process, wilcox.pvalue = extraList[[1]])

      wilcox.pvalues <- bind_rows(wilcox.pvalues, wilcoxVector) %>%
        mutate(p.adj = p.adjust(wilcox.pvalue, method = "BH"),
               significant = ifelse(p.adj <= 0.05, "*", ""))

      # extracts the chi square value for Gleason proportion differences
      chiVector <- c(Process = process, chi.pvalue = extraList[[2]])
      chi.pvalues <- bind_rows(chi.pvalues, chiVector) %>%
        mutate(p.adj = p.adjust(chi.pvalue, method = "BH"),
               significant = ifelse(p.adj <= 0.05, "*", ""))
    }

  }

  if(exists("wilcox.pvalues") && !is.null(extraList)){
    write_csv(wilcox.pvalues, here::here("TCGA_results", project, automataID, "clinical_analysis",
                                         paste0("wilcoxPSA_", project, "_", automataID, ".csv")))
  }

  if(exists("chi.pvalues") && !is.null(extraList)){
    write_csv(chi.pvalues, here::here("TCGA_results", project, automataID, "clinical_analysis",
                                      paste0("chisqGleason_", project, "_", automataID, ".csv")))
  }


  write_csv(summary, here::here("TCGA_results", project, automataID,
                                paste0("Summary_", project, "_", automataID, ".csv")))

  message("\n\n")
  message("=======================================================")
  message("\n\n")


}





