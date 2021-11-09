#' Creates job folders
#'
#' Generates a folder with the required data to run LPD and changes the working directory to it.
#' @param project Name of a TCGA project (eg. TCGA-PRAD or TCGA-BRCA).
#' @param automataID ID assigned to the analysis so the files from different analysis are not mixed.
#' @param path Filepath to the original working directory to return after finishing the function.
#' @param sigma Sigma value (function to run inside of a loop)
#' @importFrom stringr str_sub

generateJobFolder <- function(project, automataID, path, sigma, lpdStage, n_process = NULL){

  # Reassigns the working directory to the original one
  setwd(path)

  # Creates the directories for the id
  dir.create(paste0(path, "/TCGA_results/", project, "/", automataID, "/LPD_stage", lpdStage))

  # Creates a folder for each run using timestamp to avoid overwritting
  if(lpdStage == "B"){
    run_folder <- paste0(path,
                         "/TCGA_results/",
                         project,
                         "/",
                         automataID,
                         "/LPD_stage",
                         lpdStage,
                         "/run_",
                         format(Sys.time(), "%Y_%m_%d_%H_%M_%S"),
                         "_",
                         str_sub(as.numeric(Sys.time()), -3, -1),
                         "_",
                         n_process,
                         sigma,
                         "_B")
  } else{
    run_folder <- paste0(path,
                         "/TCGA_results/",
                         project,
                         "/",
                         automataID,
                         "/LPD_stage",
                         lpdStage,
                         "/run_",
                         format(Sys.time(), "%Y_%m_%d_%H_%M_%S"),
                         "_",
                         str_sub(as.numeric(Sys.time()), -3, -1),
                         "_",
                         sigma,
                         "_A")
  }

  dir.create(run_folder)
  print(run_folder)
  setwd(run_folder)

}

#' Runs the first stage of LPD
#'
#' Launches 90 jobs trying different sigma values with different number of processes.
#' @param project Name of a TCGA project (eg. TCGA-PRAD or TCGA-BRCA).
#' @param automataID ID assigned to the analysis so the files from different analysis are not mixed.
#' @param path Pathway to TCGA_results folder. Defaults to working directory.
#' @importFrom  readr read_csv
#' @importFrom dplyr select
#' @export

runLPD_A <- function(project, automataID, path = getwd()){

  # Imports the data
  data <- read_csv(paste0(path,
                          "/TCGA_results/",
                          project,
                          "/",
                          automataID,
                          "/processed_expression/topGenes_",
                          project,
                          "_",
                          automataID,
                          ".csv")) %>%
    dplyr::select(-c(1))


  # Sigmas to analyse
  sigmas <- c(-0.001, -0.005, -0.01, -0.05, -0.08,
              -0.1, -0.2, -0.4, -0.5, -0.6, -0.7, -0.8,
              -0.9, -1, -1.1, -1.2, -1.4, -1.5)

  # Looping through sigmas
  for(s in sigmas){

    # Do 5 replicates for sigma
    rep <- 1

    while(rep < 6){

      # Avoids having the same timestamp and overwriting
      Sys.sleep(0.15)
      generateJobFolder(project, automataID, path, s, lpdStage = "A")

      # Generates a config file in each folder
      sink("Config.txt")

      cat('FILENAME data.txt\n')
      cat("MINPROCESSES 2\n")
      cat("MAXPROCESSES 15\n")
      cat("ITERATIONS 1000\n")
      cat("LOAD 0\n")
      cat(paste0("RANDOMSEED ", str_sub(as.numeric(Sys.time()), -3, -1), "\n"))
      cat("MEANPRIORM 0.0\n")
      cat(paste0("SIGMAPRIOR1 "), s, "\n")
      cat("SIGMAPRIOR2 0\n")
      cat("SIGMAPRIOR3 0\n")
      cat("CROSSVALIDATE 1\n")
      cat("VERBOSE 0\n")
      cat("LOOP 1")

      sink()

      # Generates a job file in each folder
      sink("job.bsub")

      cat('#!/bin/sh\n')
      cat("#BSUB -q long-eth\n")
      cat(paste0("#BSUB -J ", "A_", automataID, "\n"))
      cat("#BSUB -oo out_LPD-%J.out\n")
      cat("#BSUB -eo err_LPD-%J.err\n\n")
      cat(paste0(". /etc/profile ",  getwd(),"\n"))
      cat("module add gsl\n")
      cat("/gpfs/afm/cancergenetics/Sergio/lpd_source/EMLPD")

      sink()

      # Saves the data as txt file and adds the row and col numbers
      txt <- sprintf("%s %s", nrow(data), ncol(data))
      tmp <- "data.txt"
      cat(txt, "\n", file = tmp)
      write.table(data, tmp, row.names = FALSE, col.names = FALSE, append = TRUE)
      unlink("tmp.txt")

      # Launches LPD
      system("bsub < job.bsub")

      rep <- rep + 1
    }
    setwd(path)
  }

}


#' Runs the first stage of LPD (ada)
#'
#' Launches 90 jobs trying different sigma values with different number of processes for ADA.
#' @param project Name of a TCGA project (eg. TCGA-PRAD or TCGA-BRCA).
#' @param automataID ID assigned to the analysis so the files from different analysis are not mixed.
#' @param path Pathway to TCGA_results folder. Defaults to working directory.
#' @importFrom  readr read_csv
#' @importFrom dplyr select
#' @export

runLPD_A_ADA <- function(project, automataID, path = getwd()){

  # Imports the data
  data <- read_csv(paste0(path,
                          "/TCGA_results/",
                          project,
                          "/",
                          automataID,
                          "/processed_expression/topGenes_",
                          project,
                          "_",
                          automataID,
                          ".csv")) %>%
    dplyr::select(-c(1))


  # Sigmas to analyse
  sigmas <- c(-0.001, -0.005, -0.01, -0.05, -0.08,
              -0.1, -0.2, -0.4, -0.5, -0.6, -0.7, -0.8,
              -0.9, -1, -1.1, -1.2, -1.4, -1.5)

  # Looping through sigmas
  for(s in sigmas){

    # Do 5 replicates for sigma
    rep <- 1

    while(rep < 6){

      # Avoids having the same timestamp and overwriting
      Sys.sleep(0.15)
      generateJobFolder(project, automataID, path, s, lpdStage = "A")

      # Generates a config file in each folder
      sink("Config.txt")

      cat('FILENAME data.txt\n')
      cat("MINPROCESSES 2\n")
      cat("MAXPROCESSES 15\n")
      cat("ITERATIONS 1000\n")
      cat("LOAD 0\n")
      cat(paste0("RANDOMSEED ", str_sub(as.numeric(Sys.time()), -3, -1), "\n"))
      cat("MEANPRIORM 0.0\n")
      cat(paste0("SIGMAPRIOR1 "), s, "\n")
      cat("SIGMAPRIOR2 0\n")
      cat("SIGMAPRIOR3 0\n")
      cat("CROSSVALIDATE 1\n")
      cat("VERBOSE 0\n")
      cat("LOOP 1")

      sink()

      # Generates a job file in each folder
      sink("job.bsub")

      cat('#!/bin/sh\n')
      cat("#SBATCH -p compute-24-96\n")
      cat("#SBATCH -t 168:00:00\n")
      cat(paste0("#SBATCH --job-name=", "A_", automataID, "\n"))
      cat("#SBATCH -o out_LPD-%J.out\n")
      cat("#SBATCH -e err_LPD-%J.err\n\n")
      cat(paste0(". /etc/profile ",  getwd(),"\n"))
      cat("module add gsl/2.2\n")
      cat("/gpfs/afm/cancergenetics/Sergio/lpd_source/EMLPD")

      sink()

      # Saves the data as txt file and adds the row and col numbers
      txt <- sprintf("%s %s", nrow(data), ncol(data))
      tmp <- "data.txt"
      cat(txt, "\n", file = tmp)
      write.table(data, tmp, row.names = FALSE, col.names = FALSE, append = TRUE)
      unlink("tmp.txt")

      # Launches LPD
      system("sbatch job.bsub")

      rep <- rep + 1
    }
    setwd(path)
  }

}




#' Estimate best number of processes and sigma value
#'
#' Creates a file with the best combinations and plots the likelihood for each one
#' @param project Name of a TCGA project (eg. TCGA-PRAD or TCGA-BRCA).
#' @param automataID ID assigned to the analysis so the files from different analysis are not mixed.
#' @param path Path of the working directory
#' @importFrom dplyr bind_rows filter
#' @importFrom plyr ddply
#' @import magrittr
#' @importFrom readr write_csv
#' @import ggplot2
#' @export

estimateParameters <- function(project, automataID, path = getwd()){

  # Arguments
  folderName <- paste0(path, "/TCGA_results/", project, "/", automataID, "/LPD_stageA/")
  sigmaValues <-  c(-0.001, -0.005, -0.01, -0.05, -0.08,
                    -0.1, -0.2, -0.4, -0.5, -0.6, -0.7, -0.8,
                    -0.9, -1, -1.1, -1.2, -1.4, -1.5)

  files_and_sigmas <- data.frame(files = c(), sigmas = c())

  # Read all the sigmas from each file
  for(sigma in sigmaValues) {
    filenames <- list.files(path = folderName, pattern = paste("(.*)", "_", sigma, "_A", sep = ""))
    temp = data.frame(files = filenames, sigmas = sigma)
    files_and_sigmas = bind_rows(files_and_sigmas, temp)
  }

  all_Likelihood_raw = data.frame(Processes = c(), Likelihood = c(), Sigma = c())

  for(filename in files_and_sigmas$files) {
    data = read.csv(paste(folderName, filename, "/LogLikelihood.txt", sep = ""), sep = "", header = F)
    colnames(data) = c("Processes", "Likelihood")
    data$Sigma = files_and_sigmas$sigmas[files_and_sigmas$files == filename]
    all_Likelihood_raw = rbind(all_Likelihood_raw, data)
  }

  summaryData = ddply(all_Likelihood_raw, c("Processes", "Sigma"), summarise, N=length(Likelihood),
                      mean=mean(Likelihood, na.rm = TRUE), sd=sd(Likelihood, na.rm = TRUE), se=(sd/sqrt(N)))
  summaryData$ci = summaryData$se * (qt((0.95/2)+0.5, summaryData$N-1))

  # Choose best combinations
  supreme <- summaryData %>%
    .[which.min(abs(.$mean)),]

  best_option <- summaryData %>%
    filter(.$mean >= (supreme$mean - supreme$sd)) %>%
    .[which.min(.$Processes),] %>%
    .[which.min(.$mean),]

  best_other_one <- summaryData %>%
    filter(Processes == best_option$Processes + 1) %>%
    filter(mean == max(mean))

  best_other_two <- summaryData %>%
    filter(Processes == best_option$Processes - 1) %>%
    filter(mean == max(mean))

  # Export it
  write_csv(bind_rows(best_option, best_other_one, best_other_two),
            paste0("TCGA_results/",
                   project, "/",
                   automataID,
                   "/lpd_parameters/bestProcessesAndSigmas_",
                   project, "_",
                   automataID, ".csv"))

  # Plot it
  pd = position_dodge(0.3)

  p = ggplot(summaryData, aes(x = as.factor(Processes), y= mean, colour = as.factor(Sigma), group = Sigma))
  q = p + geom_errorbar(aes(ymin = mean - ci, ymax = mean + ci), colour = "black", width = .1, position = pd) +
    geom_point(position = pd, size = 3) +
    ggtitle(paste(project, automataID, "Likelihood ", sep = " ")) + xlab("Processes") +
    ylab("Log-Likelihood") +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(colour='Sigma value')

  pdf(paste0("TCGA_results/", project, "/", automataID, "/lpd_parameters/Likelihood_", project, "_", automataID, ".pdf"))
  plot(q)
  dev.off()

  # Saves the plot as an object to use it in the report lately
  likelihood.plot <- q
  save(likelihood.plot, file = paste0("TCGA_results/", project, "/", automataID,
                                      "/report_plots/Likelihood_Plot_", project, "_", automataID, ".Rdata"))

  return(bind_rows(best_option, best_other_one, best_other_two))
}

#' Runs the second stage of LPD
#'
#' Launches 90 jobs trying different sigma values with different number of processes.
#'
#' @param project Name of a TCGA project (eg. TCGA-PRAD or TCGA-BRCA).
#' @param automataID ID assigned to the analysis so the files from different analysis are not mixed.
#' @param n_process Number of processes estimated as optimal
#' @param sigma Sigma value estimated as optimal.
#' @param path Pathway to TCGA_results folder. Defaults to working directory.
#' @importFrom  readr read_csv
#' @importFrom dplyr select
#' @export

launchLPD_B <- function(project, automataID, n_process, sigma, path = getwd()){

  # Imports the data
  data <- read_csv(paste0(path,
                          "/TCGA_results/",
                          project,
                          "/",
                          automataID,
                          "/processed_expression/topGenes_",
                          project,
                          "_",
                          automataID,
                          ".csv")) %>%
    dplyr::select(-c(1))


  # Only one sigma is analysed
  s <- sigma

  # Do 100 replicates for sigma
  rep <- 1

  while(rep < 101){

    # Avoids having the same timestamp and overwriting
    Sys.sleep(0.15)
    generateJobFolder(project, automataID, path, s, lpdStage = "B", n_process = n_process)

    # Generates a config file in each folder
    sink("Config.txt")

    cat('FILENAME data.txt\n')
    cat(paste0("MINPROCESSES ", n_process, "\n"))
    cat(paste0("MAXPROCESSES ", n_process, "\n"))
    cat("ITERATIONS 1000\n")
    cat("LOAD 0\n")
    cat(paste0("RANDOMSEED ", str_sub(as.numeric(Sys.time()), -3, -1), "\n"))
    cat("MEANPRIORM 0.0\n")
    cat(paste0("SIGMAPRIOR1 "), s, "\n")
    cat("SIGMAPRIOR2 0\n")
    cat("SIGMAPRIOR3 0\n")
    cat("CROSSVALIDATE 0\n")
    cat("VERBOSE 0\n")
    cat("LOOP 1")

    sink()

    # Generates a job file in each folder
    sink("job.bsub")

    cat('#!/bin/sh\n')
    cat("#BSUB -q long-eth\n")
    cat(paste0("#BSUB -J ", "B_", automataID, "\n"))
    cat("#BSUB -oo out_LPD-%J.out\n")
    cat("#BSUB -eo err_LPD-%J.err\n\n")
    cat(paste0(". /etc/profile ",  getwd(),"\n"))
    cat("module add gsl\n")
    cat("/gpfs/afm/cancergenetics/Sergio/lpd_source/EMLPD")

    sink()

    # Saves the data as txt file and adds the row and col numbers
    txt <- sprintf("%s %s", nrow(data), ncol(data))
    tmp <- "data.txt"
    cat(txt, "\n", file = tmp)
    write.table(data, tmp, row.names = FALSE, col.names = FALSE, append = TRUE)
    unlink("tmp.txt")

    # Launches LPD
    system("bsub < job.bsub")

    rep <- rep + 1
  }
  setwd(path)
}


#' Runs the second stage of LPD (AADA)
#'
#' Launches 90 jobs trying different sigma values with different number of processes. COmpatible with ADA.
#'
#' @param project Name of a TCGA project (eg. TCGA-PRAD or TCGA-BRCA).
#' @param automataID ID assigned to the analysis so the files from different analysis are not mixed.
#' @param n_process Number of processes estimated as optimal
#' @param sigma Sigma value estimated as optimal.
#' @param path Pathway to TCGA_results folder. Defaults to working directory.
#' @importFrom  readr read_csv
#' @importFrom dplyr select
#' @export

launchLPD_B_ADA <- function(project, automataID, n_process, sigma, path = getwd()){

  # Imports the data
  data <- read_csv(paste0(path,
                          "/TCGA_results/",
                          project,
                          "/",
                          automataID,
                          "/processed_expression/topGenes_",
                          project,
                          "_",
                          automataID,
                          ".csv")) %>%
    dplyr::select(-c(1))


  # Only one sigma is analysed
  s <- sigma

  # Do 100 replicates for sigma
  rep <- 1

  while(rep < 101){

    # Avoids having the same timestamp and overwriting
    Sys.sleep(0.15)
    generateJobFolder(project, automataID, path, s, lpdStage = "B", n_process = n_process)

    # Generates a config file in each folder
    sink("Config.txt")

    cat('FILENAME data.txt\n')
    cat(paste0("MINPROCESSES ", n_process, "\n"))
    cat(paste0("MAXPROCESSES ", n_process, "\n"))
    cat("ITERATIONS 1000\n")
    cat("LOAD 0\n")
    cat(paste0("RANDOMSEED ", str_sub(as.numeric(Sys.time()), -3, -1), "\n"))
    cat("MEANPRIORM 0.0\n")
    cat(paste0("SIGMAPRIOR1 "), s, "\n")
    cat("SIGMAPRIOR2 0\n")
    cat("SIGMAPRIOR3 0\n")
    cat("CROSSVALIDATE 0\n")
    cat("VERBOSE 0\n")
    cat("LOOP 1")

    sink()

    # Generates a job file in each folder
    sink("job.bsub")

    cat('#!/bin/sh\n')
    cat("#SBATCH -p compute-24-96\n")
    cat("#SBATCH -t 168:00:00\n")
    cat(paste0("#SBATCH --job-name=", "B_", automataID, "\n"))
    cat("#SBATCH -o out_LPD-%J.out\n")
    cat("#SBATCH -e err_LPD-%J.err\n\n")
    cat(paste0(". /etc/profile ",  getwd(),"\n"))
    cat("module add gsl/2.2\n")
    cat("/gpfs/afm/cancergenetics/Sergio/lpd_source/EMLPD")
    sink()

    # Saves the data as txt file and adds the row and col numbers
    txt <- sprintf("%s %s", nrow(data), ncol(data))
    tmp <- "data.txt"
    cat(txt, "\n", file = tmp)
    write.table(data, tmp, row.names = FALSE, col.names = FALSE, append = TRUE)
    unlink("tmp.txt")

    # Launches LPD
    system("sbatch job.bsub")

    rep <- rep + 1
  }
  setwd(path)
}

#' Launches jobs with best combinations for LPD-B
#'
#' Launches 90 jobs trying different sigma values with different number of processes.
#' @param project Name of a TCGA project (eg. TCGA-PRAD or TCGA-BRCA).
#' @param automataID ID assigned to the analysis so the files from different analysis are not mixed.
#' @param n_process Number of processes estimated as optimal
#' @param sigma Sigma value estimated as optimal.
#' @param path Pathway to TCGA_results folder. Defaults to working directory.
#' @importFrom  readr read_csv
#' @export

runLPD_B <- function(project, automataID, path = getwd()){

  data <- read_csv(paste0(path, "/TCGA_results/", project, "/", automataID,
                          "/lpd_parameters/bestProcessesAndSigmas_", project, "_", automataID, ".csv"))

  for(i in 1:3){
    n_jobs <- length(system(
      paste0("bjobs | grep ", substr(automataID, nchar(automataID) - 9, nchar(automataID))),
      intern = TRUE))

    while(n_jobs > 130){
      message("Not enough available jobs. Trying again in 30 minutes...")
      Sys.sleep(1300)
      while(length(system("bjobs", intern = TRUE)) > 145){
        message("Not enough available jobs. Trying again in 30 minutes...")
        Sys.sleep(1300)
      }
      n_jobs <- length(system(
        paste0("bjobs | grep ", substr(automataID, nchar(automataID) - 9, nchar(automataID))),
        intern = TRUE))
    }

    launchLPD_B(project = project, automataID = automataID, n_process = data$Processes[i],
                sigma = data$Sigma[i], path = path)
  }

}



#' Launches jobs with best combinations for LPD-B (ada)
#'
#' Launches 90 jobs trying different sigma values with different number of processes.
#' @param project Name of a TCGA project (eg. TCGA-PRAD or TCGA-BRCA).
#' @param automataID ID assigned to the analysis so the files from different analysis are not mixed.
#' @param n_process Number of processes estimated as optimal
#' @param sigma Sigma value estimated as optimal.
#' @param path Pathway to TCGA_results folder. Defaults to working directory.
#' @importFrom  readr read_csv
#' @export
runLPD_B_ADA <- function(project, automataID, path = getwd()){

  data <- read_csv(paste0(path, "/TCGA_results/", project, "/", automataID,
                          "/lpd_parameters/bestProcessesAndSigmas_", project, "_", automataID, ".csv"))

  for(i in 1:3){
    n_jobs <- length(system(
      paste0("squeue -u ksa18wyu | grep ", substr(automataID, 1, 6)),
      intern = TRUE))

    while(n_jobs > 130){
      message("Not enough available jobs. Trying again in 45 minutes...")
      Sys.sleep(2700)
      while(length(system("squeue - u ksa18wyu", intern = TRUE)) > 145){
        message("Not enough available jobs. Trying again in 45 minutes...")
        Sys.sleep(2700)
      }
      n_jobs <- length(system(
        paste0("squeue -u ksa18wyu | grep ", substr(automataID, 1, 6)),
        intern = TRUE))
    }

    launchLPD_B_ADA(project = project, automataID = automataID, n_process = data$Processes[i],
                sigma = data$Sigma[i], path = path)
  }

}
