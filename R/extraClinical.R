#' Extra clinical analysis
#'
#' Performs further analysis for some specific cancers, so far it's only available
#' for TCGA-PRAD.
#'
#' @param project .
#' @param automataID .
#' @param clinicalData Import the clinical data from the main clinical analysis function
#' @param process LPD process
#' @importFrom dplyr mutate filter select group_by tally
#' @importFrom survival coxph Surv
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom tidyr drop_na spread
#' @export

extraClinical <- function(project, automataID, clinicalData, process){

  if(project != "TCGA-PRAD"){
    return(NULL)
  }

  # Creayes the dataframe
  clinicalDF <- data.frame(Sample_ID = clinicalData$bcr_patient_barcode,
                           LPD_Group = clinicalData$Max_LPD,
                           Time = clinicalData$days_to_first_biochemical_recurrence,
                           FollowUp = clinicalData$days_to_last_followup,
                           Event = clinicalData$biochemical_recurrence,
                           SelectedGroup = clinicalData$SelectedGroup,
                           Gleason = clinicalData$stage_event_gleason_grading,
                           PSA = as.numeric (clinicalData$stage_event_psa)) %>%
    mutate(Gleason = substr(Gleason, 0, 1),
           Gleason = ifelse(Gleason == 1, 10, Gleason),
           Gleason = ifelse(Gleason < 8, "Low", 'High')) %>%
    mutate(Time = ifelse(is.na(Time), FollowUp, Time)) %>%
    dplyr::filter(is.na(Event) == FALSE, is.na(Time) == FALSE) %>%
    dplyr::select(-FollowUp) %>%
    mutate(Event = ifelse(Event == "YES", TRUE, FALSE))

  # Runs cox
  coxTable <- coxph(Surv(Time, Event) ~ SelectedGroup+PSA+Gleason, data = clinicalDF)

  coxOutcome <- summary(coxTable)$coefficients %>%
    as.data.frame() %>%
    rownames_to_column("Variable")

  write_csv(coxOutcome, here::here("TCGA_results", project, automataID, "clinical_analysis",
                                   paste0("coxAnalysis_", project, "_", process, "_", automataID, ".csv")))

  coxOutcome

  # Mann-whitney for PSA
  mann.matrix <- clinicalDF %>%
    select(SelectedGroup, PSA) %>%
    drop_na()

  wilcoxOutcome <- wilcox.test(PSA ~ SelectedGroup, data = mann.matrix)

  # Cotingency table for Gleason
  contingecyGleason <- clinicalDF %>%
    group_by(SelectedGroup, Gleason) %>%
    tally() %>%
    spread(key = Gleason, value = n) %>%
    column_to_rownames("SelectedGroup") %>%
    as.matrix()

  write.csv(contingecyGleason, here::here("TCGA_results", project, automataID, "clinical_analysis",
                                   paste0("contingecyGleason_", project, "_", process, "_", automataID, ".csv")))

  chisqOutcome <- chisq.test(contingecyGleason)


  outcome <- list(Wilcox = wilcoxOutcome$p.value, ChiSq = chisqOutcome$p.value)

  outcome



}
