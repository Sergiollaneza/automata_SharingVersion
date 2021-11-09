#' Sets the same order of gamma values columns in all the runs
#'
#' Takes the first run of a combination and sorts the rest of the runs according to that one. In order
#' to identify the corresponding column, Spearman's correlation is used.
#'
#' @param project Name of the project to analyse according to the TCGA (e.g. TCGA-PRAD)
#' @param automataID Custom ID tag
#' @import here
#' @importFrom  readr read_csv write_csv
#' @importFrom dplyr bind_rows filter group_by arrange ungroup
#' @importFrom plyr adply

sortGamma <- function(project, automataID){

  # Read the best three n-process + sigma
  processAndSigma <- read_csv(here::here("TCGA_results", project, automataID, "lpd_parameters",
                                         paste0("bestProcessesAndSigmas_", project, "_",
                                                automataID, ".csv")))

  for(combination in 1:nrow(processAndSigma)){

    n_process <- processAndSigma$Processes[combination]
    sigma <- processAndSigma$Sigma[combination]

    # I need to import all the gamma values
    subfolders <- list.dirs(here::here("TCGA_results", project, automataID, "LPD_stageB")) %>%
      .[grep(paste0("_", n_process, sigma), .)]

    gammaList <- list()
    count <- 1

    for(folder in subfolders){

      if(file.exists(paste0(folder, "/Gamma.txt"))){

        gammaList[[count]] <- read.csv(paste0(folder, "/Gamma.txt"), sep = "", header = F)

        count <- 1 + count

      } else{
        next()
      }
    }


    # Let's start with Dan super code


    # Itinerates through all the runs
    for(run in 1:length(gammaList)){

      message("oo Analysing run:", subfolders[[run]])

      correlationResults <- tibble()

      # Creates a dataframe with the Rho estimation of each combination of columns
      for(i in 1:ncol(gammaList[[1]])){
        for(j in 1:ncol(gammaList[[run]])){
          rhoEstimate <- cor.test(gammaList[[1]][, i], gammaList[[run]][, j], method = "spearman")$estimate
          tempVector <- c(FirstRunColumn = i, OtherRunColumn = j, runNumber = run, Rho = rhoEstimate)
          correlationResults <- bind_rows(correlationResults, tempVector)
        }
      }

      # Filter the correlations lower than 0.1
      correlationResFiltered <- dplyr::filter(correlationResults, Rho.rho > 0.1)

      # Checks that all the columns from the first run have at least one rho that has
      # passed the threshold
      for(k in base::setdiff(1:ncol(gammaList[[1]]), correlationResFiltered$FirstRunColumn)){

        # If they don't, it puts back all of them even if they don't pass the threshold
        correlationResFiltered <- bind_rows(correlationResFiltered,
                                            dplyr::filter(correlationResults, FirstRunColumn == k))
      }

      # Each first-run column is separated into a different list to calculate permutations
      # We only keep the columns from the other run (so to which one did the first one matched)
      # It generates a table with all possible combinations of columns that still have a good match
      permutations <- split(correlationResFiltered[1:2],
                            f = correlationResFiltered$FirstRunColumn) %>%
        lapply(., function(df){
          df$OtherRunColumn
        }) %>%
        expand.grid()

      # Some permutations repeat columns, so they have to be removed
      index <- adply(permutations, 1, function(row){
        data.frame(include = length(setdiff(1:ncol(gammaList[[1]]), row)) == 0)
      })

      permutationsFiltered <- permutations[index$include, ]

      # Sometimes there are not permutation without repetting
      if(nrow(permutationsFiltered) > 0){

        # Loop to compare which permutations gives the best results overall
        metricMax <- 0

        for(rowNumber in 1:nrow(permutationsFiltered)){
          reference <- gammaList[[1]]
          testRun <- gammaList[[run]][,as.integer(permutationsFiltered[rowNumber,])]

          # Loop to calculate the pairwise correlation between columns of run 1 and test run
          corResult <- 0
          for(i in 1:ncol(gammaList[[1]])){
            tempRho <- cor.test(reference[,i], testRun[,i], method = "spearman")$estimate
            corResult <- corResult + tempRho
          }

          # The biggest correlation is saved
          if(corResult > metricMax){
            selectedPermutation <-permutationsFiltered[rowNumber,]
            metricMax <- corResult
          }

          # Reorder the gamma matrix
          gammaList[[run]] <- gammaList[[run]][,as.integer(selectedPermutation)]

        }
      } else{
        # When there is not a permutation without repeating, then each column is assigned to
        # the one that has the highest correlation

        bestMatch <- tibble(FirstRunColumn = numeric(),
                            OtherRunColumn = numeric(),
                            runNumber = numeric(),
                            Rho.rho = numeric())

        while(nrow(bestMatch) != ncol(gammaList[[1]])){

          tempMatch <- correlationResults %>%
            dplyr::filter(!(FirstRunColumn %in% bestMatch$FirstRunColumn  ),
                          !(OtherRunColumn %in% bestMatch$OtherRunColumn)) %>%
            group_by(FirstRunColumn) %>%
            arrange(desc(Rho.rho)) %>%
            slice(1) %>%
            ungroup() %>%
            group_by(OtherRunColumn) %>%
            arrange(desc(Rho.rho)) %>%
            slice(1)

          bestMatch <- bind_rows(bestMatch, tempMatch) %>%
            arrange(FirstRunColumn)
        }


        gammaList[[run]] <- gammaList[[run]][,as.integer(bestMatch$OtherRunColumn)]


      }

      write.csv(gammaList[[run]],
                paste0(subfolders[run], "/GammaSorted.txt"))


    }
  }

}



#' Compare gammas to select the best run of each combination
#'
#' Checks each one of the 100 gamma files and calculates a mean gamma. The file with the lowest distance
#' between the mean and its values is chosen as the most representative.
#'
#' @param project Name of the project to analyse according to the TCGA (e.g. TCGA-PRAD)
#' @param automataID Custom ID tag
#' @import here
#' @importFrom  tibble tibble
#' @importFrom  readr read_csv write_csv
#' @importFrom dplyr bind_rows
#' @import magrittr
#' @importFrom abind abind

compareGammas <- function(project, automataID){

  representants <- tibble()


  # Read the best three n-process + sigma
  processAndSigma <- read_csv(here::here("TCGA_results", project, automataID, "lpd_parameters",
                                         paste0("bestProcessesAndSigmas_", project, "_",
                                                automataID, ".csv")))

  for(row in 1:nrow(processAndSigma)){

    n_process <- processAndSigma$Processes[row]
    sigma <- processAndSigma$Sigma[row]

    folderList <- list.dirs(here::here("TCGA_results", project, automataID, "LPD_stageB")) %>%
      .[grep(paste0("_", n_process, sigma), .)]

    index = 1
    gammaList <- list()

    for(subFolder in folderList){

      # Read gamma of the subfolder

      if(file.exists(paste0(subFolder, "/GammaSorted.txt"))){
        gammaData <- read.csv(paste0(subFolder, "/GammaSorted.txt"),  header = T, row.names = 1)
      } else{
        next()
      }

      nameVector <- c(paste("LPD_",(1:length(gammaData)),sep=""))
      colnames(gammaData) <- nameVector

      gammaList[[index]] <- gammaData

      index <-index + 1

    }

    meanGammas <- apply(abind(gammaList, along = 3), 1:2, mean)

    diffGammas <- tibble()
    # Calculate which one has the lowest total difference to the mean
    for(subFolder in folderList){

      # Read gamma of the subfolder
      if(file.exists(paste0(subFolder, "/GammaSorted.txt"))){
        gammaData <- read.csv(paste0(subFolder, "/GammaSorted.txt"),  header = T, row.names = 1)
      } else{
        next()
      }

      diff <- abs(meanGammas - gammaData)
      totalDiff <- sum(diff)

      tempVector <- c(ID = subFolder, Distance = totalDiff, Process = n_process, Sigma = sigma)

      diffGammas <- bind_rows(diffGammas, tempVector) %>%
        arrange(Distance)
    }

    selected_run <- diffGammas[1,]

    representants <- bind_rows(representants, selected_run[1,])

  }

  write_csv(representants, here::here("TCGA_results", project, automataID, "lpd_parameters",
                                      paste0("bestRuns_", project, "_", automataID, ".csv")))

}

#' Calculates correlation and chooses best combination
#'
#' Calculates the correlation for each combination and selects the one with closest Pearson to 0 and
#' p-value lower than 0.05 as the best.
#' @param project Name of the project to analyse according to the TCGA (e.g. TCGA-PRAD)
#' @param automataID Custom ID tag
#'
#' @importFrom readr read_csv write_csv read_delim
#' @import here
#' @importFrom dplyr tibble mutate bind_rows select_if group_by summarize filter top_n select
#' @import magrittr
#'

calculateCorrelation <- function(project, automataID, genes){

  correlations <- tibble()

  # Read the best run for each combination
  bestRuns <- read_csv(here::here("TCGA_results", project, automataID, "lpd_parameters",
                                  paste0("bestRuns_", project, "_", automataID, ".csv")))


  for(row in 1:nrow(bestRuns)){
    folder <- bestRuns$ID[row]
    n_process <- bestRuns$Process[row]
    sigma <- bestRuns$Sigma[row]

    # Read the Means.txt file of the selected run and transpose it
    meanData <- read_delim(paste0(folder, "/Means.txt"), delim = " ", col_names = FALSE) %>%
      select_if(~sum(!is.na(.)) > 0) %>%
      t()

    colnames(meanData) <- colnames(genes)[-1]

    # Calculate mean and sd for each column
    means <- apply(meanData, 2 , mean)
    sds <- apply(meanData, 2 , sd)

    # Calculate how different is each value for its mean column value and divides by SD
    meanNormal <- meanData %>%
      sweep(2, means) %>%
      sweep(2, sds, "/")

    # Calculates correlation between each cluster (row) with the other
    correlation_temp <- tibble()

    for(reference in 1:nrow(meanNormal)){
      for(target in 1:nrow(meanNormal)){
        if(reference != target && reference < target){
          corrValue <- cor(meanNormal[reference,], meanNormal[target,], method = "pearson")
          corrTest <- cor.test(meanNormal[reference,], meanNormal[target,], method = "pearson")$p.value
          temp_df <- tibble(Reference = meanNormal[reference,],
                            Reference_process = reference,
                            Target = meanNormal[target,],
                            Target_process = target,
                            Pearson = corrValue,
                            pVal = corrTest)
          correlation_temp <- bind_rows(correlation_temp, temp_df)
        }
      }
    }

    # Summarise the correlations
    correlations_summ <- correlation_temp %>%
      dplyr::group_by(Reference_process, Target_process) %>%
      summarize(Mean_reference = mean(Reference),
                Mean_target = mean(Target),
                Pearson = mean(Pearson),
                pVal = mean(pVal),
                n_process = n_process,
                sigma = sigma) %>%
      mutate(Folder = folder)


    correlations <- bind_rows(correlations, correlations_summ)
    write_csv(correlations, here::here("TCGA_results", project, automataID, "lpd_parameters",
                                       paste0("correlations_", project, "_", automataID, ".csv")))

  }

  bestCombination <- correlations %>%
    filter(pVal < 0.05) %>%
    mutate(Pearson.square = (Pearson ^ 2)) %>%
    group_by(n_process, sigma, Folder) %>%
    summarise(Mean.square.Pearson = mean(Pearson.square)) %>%
    ungroup()  %>%
    top_n(-1, Mean.square.Pearson)

  write_csv(bestCombination, here::here("TCGA_results", project, automataID, "lpd_parameters",
                                        paste0("bestCombination_", project, "_", automataID, ".csv")))

}


#' Assigns each sample to an LPD group according to the gamme value
#'
#' Assigns each sample to an LPD group according to the gamme value
#' @param project Name of the project to analyse according to the TCGA (e.g. TCGA-PRAD)
#' @param automataID Custom ID tag
#' @importFrom readr read_csv write_csv
#' @import here
#' @importFrom dplyr mutate arrange select
#' @importFrom reshape2 melt
#' @import ggplot2
#' @import magrittr
#'

assignGammas <- function(project, automataID, genes){

  bestCombo <- read_csv(here::here("TCGA_results", project, automataID, "lpd_parameters",
                                   paste0("bestCombination_", project, "_", automataID, ".csv")))

  folder <- bestCombo$Folder
  n_process <- bestCombo$n_process
  sigma <- bestCombo$sigma

  gammaData <- read.csv(paste0(folder, "/Gamma.txt"), sep = "", header = F)
  nameVector <- c(paste("LPD_",(1:length(gammaData)),sep=""))
  colnames(gammaData) <- nameVector

  # Normalise the values
  gamma_values <- gammaData / rowSums(gammaData)
  # Adds sample name
  gamma_values$Aliquot = gsub("\\.", "-", genes$X1)

  # Fix the factor levels to keep the order in the plot
  signalRows <- grep("Signal_", gamma_values$Aliquot)
  gamma_values$Aliquot[signalRows] <- substr(gamma_values$Aliquot[signalRows], 8, 35)

  gamma_values$Max_LPD <- colnames(gamma_values)[max.col(gamma_values[,1:n_process], ties.method="first")]
  gamma_values <- arrange(gamma_values, Max_LPD)
  gamma_values$Index <- 1:nrow(gamma_values)
  gamma_values$Sample <- substr(gamma_values$Aliquot, 1, 16)
  gamma_values$Aliquot <- factor(gamma_values$Aliquot, levels = gamma_values$Aliquot)

  write_csv(gamma_values,here::here("TCGA_results", project, automataID, "lpd_parameters",
                       paste0("GammaValues_FOR4_", project, "_", automataID, ".csv")))

  write_csv(dplyr::select(gamma_values, Sample, Aliquot, Max_LPD),
            here::here("TCGA_results", project, automataID, "lpd_parameters",
                       paste0("Sample_assigned_FOR4_", project, "_", automataID, ".csv")))
  # Reshapes to long format
  # Aliquot will replace the sample column to avoid overlapping
  gamma_melt <- gamma_values %>%
    select(-c(Sample)) %>%
    melt(id.vars = c("Aliquot", "Max_LPD", "Index"))

  colnames(gamma_melt) <- c("Sample", "Max_LPD", "Index", "LPD_Process", "Gamma")


  p <- ggplot(gamma_melt, aes(x = reorder(Sample, Index), y = Gamma, fill = Max_LPD)) +
    geom_bar(stat = "identity") +
    facet_grid(LPD_Process ~ .) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text.x = element_blank()) +
    xlab("Sample") +
    ggtitle(paste("Gamma plot for", project, automataID, n_process, sigma, sep = " ")) +
    scale_fill_brewer(palette = ifelse(length(unique(gamma_melt$Max_LPD)) < 9, "Set2", "Paired"))

  pdf(here::here("TCGA_results", project, automataID, "lpd_parameters",
                 paste0("GammaPlot_FOR4_", project, "_", n_process, "-", sigma, "_", automataID, ".pdf")))
  plot(p)
  dev.off()


  gammaSample.plot <- p
  save(gammaSample.plot, file = paste0("TCGA_results/", project, "/", automataID,
                                       "/report_plots/gammaSample_Plot_", project, "_", automataID, ".Rdata"))

}


#' Postprocessing
#'
#' Process each run and select the best one across the three combinationss and assigns each sample
#' to the LPD group to which they are more represented.
#' @param project Name of the project to analyse according to the TCGA (e.g. TCGA-PRAD)
#' @param automataID Custom ID tag
#'
#' @importFrom readr read_csv
#' @import here
#' @export

postprocess <- function(project, automataID){

  message("o Postprocess...")

  # Reads the top 500 genes
  topGenes <- read_csv(here::here("TCGA_results", project, automataID, "processed_expression",
                                  paste0("topGenes_", project, "_", automataID, ".csv")))

  message("oo Sorting gammas...")
  sortGamma(project = project, automataID = automataID)
  message("oo Comparing gammas...")
  compareGammas(project = project, automataID = automataID)
  message("oo Calculating correlations...")
  calculateCorrelation(project = project, automataID = automataID, genes = topGenes)
  message("oo Assigning gammas...")
  assignGammas(project = project, automataID = automataID, genes = topGenes)
}
