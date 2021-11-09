createOverviewRMD <- function(project, automataID, workingPath){
  

# DATA EXTRACTED FROM BARCODES --------------------------------------------
  sink(paste0(workingPath, "/TCGA_results/", project, "/", automataID,
              "/Report_RMD/Childs/Barcodes/overview_barcodes.Rmd"))
  
  cat('\n```{r functions_overview_barcodes, include=FALSE}\n')
  cat('n_samples <- nrow(params$barcodes)\n')
  cat('n_patients <- length(unique(params$barcodes$Participant))\n')
  cat('```\n\n\n')
  cat('Data obtained from Barcodes: TITLE {data-height=50}\n')
  cat('-----------------------------------------------------------------------\n\n')
  cat('#### Analysis details \n\n')
  cat('Data obtained from Barcodes: INFO\n')
  cat('-----------------------------------------------------------------------\n\n')
  cat('### Project type and automata ID\n\n')
  cat('```{r}\n')
  cat('valueBox(value = paste0("Project: ", params$project), caption = paste0("Automata ID: ",
      params$automataID), icon = "fa-info-circle", color = "#2EA587",  href = paste0( params$path, 
      "/TCGA_results/", params$project, "/", params$automataID, 
      "/Report/Automata_Report_01.html#barcodes"))\n')
  cat('```\n\n')
  cat('### Samples available\n\n')
  cat('```{r}\n')
  cat('valueBox(value = paste0(n_samples, " samples available"), caption = paste0("for ", params$method, 
      " from ", n_patients, " patients"), icon = "fa-user-friends", color = "#2EA587", 
      href = paste0( params$path, "/TCGA_results/", params$project, "/", params$automataID, 
      "/Report/Automata_Report_01.html#barcodes"))\n')
  cat("```\n\n\n")
  cat('-----------------------------------------------------------------------\n')
  cat("\n")
  sink()
  
  
  # DATA EXTRACTED FROM BARCODES --------------------------------------------
  sink(paste0(workingPath, "/TCGA_results/", project, "/", automataID,
              "/Report_RMD/Childs/LPD_parameters/overview_LPD.Rmd"))
  
  # cat('Data obtained from LPD: TITLE {data-height=50}\n')
  # cat('-----------------------------------------------------------------------\n\n')
  # cat('#### Test \n\n')
  # cat('Data obtained from LPD: INFO\n')
  # cat('-----------------------------------------------------------------------\n\n')
  cat('### Number of processes and sigma \n\n')
  cat('```{r}\n')
  cat('valueBox(value = paste0("Nº processes: ", params$bestCombination$n_process), caption = paste0("Sigma value: ",
      params$bestCombination$sigma), icon = "fa-chart-pie", color = "#2EA587",  href = paste0( params$path, 
      "/TCGA_results/", params$project, "/", params$automataID, 
      "/Report/Automata_Report_01.html#lpd-parameters"))\n')
  cat("```\n\n\n")
  cat('-----------------------------------------------------------------------\n')
  cat("\n")
  sink()
  
  
  
  
}
