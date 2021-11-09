#' Runs a dendrogram
#'
#' Runs a dendogram to compare its results with the LPD clustering using the top 500 genes per sample
#' @param project Name of a TCGA project (eg. TCGA-PRAD or TCGA-BRCA).
#' @param automataID ID assigned to the analysis so the files from different analysis are not mixed.
#' @importFrom readr read_csv
#' @importFrom dplyr mutate rename left_join
#' @import magrittr
#' @import dendextend
#' @import RColorBrewer
#' @importFrom tibble column_to_rownames rownames_to_column
#' @importFrom tidyr pivot_wider
#' @export



runDendrogram <- function(project, automataID){

  # Dependencies ------------

  # Read the top genes that are used as input for LPD
  topGenes <- read_csv(here::here("TCGA_results", project, automataID, "processed_expression",
                                  sprintf("topGenes_%s_%s.csv", project, automataID))) %>%
    column_to_rownames("X1") %>%
    as.matrix()

  # Reads the gamma asssigned samples from LPD
  gammaValues <- read_csv(here::here("TCGA_results", project, automataID, "lpd_parameters",
                                     paste0("Sample_assigned_", project, "_", automataID, ".csv")))

  # Extracts all the processes found for the dataset
  processes <- unique(gammaValues$Max_LPD)


  # Processing -------------

  # Calculates the distance across samples
  distance <- dist(topGenes)
  hclusters <- hclust(distance, method = "complete")

  # Assigns each sample to a cluster
  assignedSamples <- cutree(hclusters, k = length(processes)) %>%
    as.data.frame() %>%
    rownames_to_column("Aliquot") %>%
    rename(hClust = ".") %>%
    mutate(Aliquot = gsub("\\.", "-", Aliquot)) %>%
    left_join(gammaValues)

  write_csv(assignedSamples, here::here("TCGA_results", project, automataID, "dendrogram",
                                        paste0("dendrogram_", project, "_", automataID, ".csv")))

  # Creates a table of proportions of hClusters per LPDclusters
  assignedProportions <- assignedSamples %>%
    group_by(Max_LPD, hClust) %>%
    tally() %>%
    pivot_wider(names_from = "hClust", values_from = "n") %>%
    replace(is.na(.), 0) %>%
    rowwise() %>%
    mutate(TotalSamples = sum(c_across("1":"2")),
           across(where(is.numeric), ~round(.x/TotalSamples, 4)))

  write_csv(assignedProportions, here::here("TCGA_results", project, automataID, "dendrogram",
                                            paste0("subtypeComparisonProportionLPD_", project, "_", automataID, ".csv")))

  # Creates the opposite, a table of proportions of LPDclusters per hClusters
  assignedProportions2 <- assignedSamples %>%
    group_by(hClust, Max_LPD) %>%
    tally() %>%
    pivot_wider(names_from = "Max_LPD", values_from = "n") %>%
    mutate(hClust = as.character(hClust)) %>%
    replace(is.na(.), 0) %>%
    rowwise() %>%
    mutate(TotalSamples = sum(c_across("LPD_1":max(processes))),
           across(where(is.numeric), ~round(.x/TotalSamples, 4)))

  write_csv(assignedProportions2, here::here("TCGA_results", project, automataID, "dendrogram",
                                             paste0("subtypeComparisonProportionhClusters_", project, "_", automataID, ".csv")))


  # Plots the dendogram --------------

  # Creates the dendrogram using Rcolorbrewer colours
  dendro <- hclusters %>%
    as.dendrogram() %>%
    set("branches_k_color", brewer.pal(length(processes),"Set1"), k = length(processes)) %>%
    set("labels_col", "white")

  # Creates a bar matching the colours of the dendrogram
  H.Cluster <- cutree(dendro, k = length(processes), order_clusters_as_data = FALSE) %>%
    as_tibble() %>%
    mutate(Colour = case_when(value == 1 ~ brewer.pal(9,"Set1")[1],
                              value == 2 ~ brewer.pal(9,"Set1")[2],
                              value == 3 ~ brewer.pal(9,"Set1")[3],
                              value == 4 ~ brewer.pal(9,"Set1")[4],
                              value == 5 ~ brewer.pal(9,"Set1")[5],
                              value == 6 ~ brewer.pal(9,"Set1")[6],
                              value == 7 ~ brewer.pal(9,"Set1")[7],
                              value == 8 ~ brewer.pal(9,"Set1")[8],
                              value == 9 ~ brewer.pal(9,"Set1")[9],
                              TRUE ~ "black")) %>%
    dplyr::pull(Colour)

  # Creates a bar with the colours of the LPD clusters
  LPD.Group <- assignedSamples[order.dendrogram(dendro),] %>%
    mutate(Colour = case_when(Max_LPD == "LPD_1" ~ brewer.pal(9,"Set3")[1],
                              Max_LPD == "LPD_2" ~ brewer.pal(9,"Set3")[2],
                              Max_LPD == "LPD_3" ~ brewer.pal(9,"Set3")[3],
                              Max_LPD == "LPD_4" ~ brewer.pal(9,"Set3")[4],
                              Max_LPD == "LPD_5" ~ brewer.pal(9,"Set3")[5],
                              Max_LPD == "LPD_6" ~ brewer.pal(9,"Set3")[6],
                              Max_LPD == "LPD_7" ~ brewer.pal(9,"Set3")[7],
                              Max_LPD == "LPD_8" ~ brewer.pal(9,"Set3")[8],
                              Max_LPD == "LPD_9" ~ brewer.pal(9,"Set3")[9],
                              TRUE ~ brewer.pal(10,"Set3")[10])) %>%
    dplyr::pull(Colour)

  # Binds them
  sideBars <- cbind(LPD.Group, H.Cluster)

  # Plots the dendrogram and the bars + a legend of the LPD groups
  pdf(here::here("TCGA_results", project, automataID, "dendrogram",
                 paste0("dendrogram_", project, "_",  automataID, ".pdf")),
      width = 16, height = 11)
  par(mar = c(10,5,2,2))
  plot(dendro)
  colored_bars(colors = sideBars)
  legend("topleft", legend = levels(as.factor(assignedSamples$Max_LPD[order.dendrogram(dendro)])),
         fill = c(brewer.pal(9,"Set3")[1], brewer.pal(9,"Set3")[2], brewer.pal(9,"Set3")[3],
                  brewer.pal(9,"Set3")[4], brewer.pal(9,"Set3")[5], brewer.pal(9,"Set3")[6],
                  brewer.pal(9,"Set3")[7], brewer.pal(9,"Set3")[8], brewer.pal(9,"Set3")[9])[1:length(processes)])



  dev.off()


}



