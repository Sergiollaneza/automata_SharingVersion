generatePathways <- function(project, automataID, process, workingPath){
  
  sink(paste0(workingPath, "/TCGA_results/", project, "/", automataID,
              "/Report_RMD/Childs/Signatures/page_pathways_", process, ".Rmd"))
  
  # Creates title
  cat('Signature', process, '{data-navmenu="Pathway analysis"}\n')
  cat('=================================================================\n\n')
  
  
  # CREATES CHILD FOR DIFF EXPRESSED GENES ---------------------------
    if(file.exists(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/expression_analysis/pathway_analysis/",
                          "egoAnalysis_overexpressedGenes_LPD_", process, "_", automataID, ".csv"))){
      
      cat('Data obtained from oVEREXPRESSED GENES: INFO--pathways\n')
      cat('-----------------------------------------------------------------------\n\n')
      
      cat("### { .mainText}\n\n")
      cat("#### Overexpressed genes \n\n")
      
      cat("Biological pathways affected by overexpressed genes for signature", process, "according to the differential gene expression analysis.\n\n")

      cat('Data obtained from DEG: PATHWAY ANALYSIS OVEREXPRESSED {.tabset} \n')
      cat('-----------------------------------------------------------------------  \n\n')

      cat("### Gene ontology pathways \n\n")
      cat("```{r table_GO_Over_LPD_", process, ", echo=FALSE}\n", sep = "")
      cat("DT::datatable(params$DEgoOver_List[[", process, "]], rownames = FALSE, options = list(columnDefs = list(list(className = 'dt-center', targets = '_all'))))\n")
      cat("```\n")
      cat("\n\n")

      cat("### KEGG pathways\n\n")
      cat("```{r table_KEGG_Over_LPD_", process, ", echo=FALSE}\n", sep = "")
      cat("DT::datatable(params$DEkeggOver_List[[", process, "]], rownames = FALSE, options = list(columnDefs = list(list(className = 'dt-center', targets = '_all'))))\n")
      cat("```\n")
      cat("\n\n")

      cat("### Cancer hallmarks \n\n")
      cat("```{r table_hallmarks_Over_LPD_", process, ", echo=FALSE}\n", sep = "")
      cat("DT::datatable(params$DEHallmarkOver_List[[", process, "]], rownames = FALSE, options = list(columnDefs = list(list(className = 'dt-center', targets = '_all'))))\n")
      cat("```\n")
      cat("\n\n")

    }

    if(file.exists(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/expression_analysis/pathway_analysis/",
                          "egoAnalysis_underexpressedGenes_LPD_", process, "_", automataID, ".csv"))){
      
      
      cat('Data obtained from underEXPRESSED GENES: INFO--pathways\n')
      cat('-----------------------------------------------------------------------\n\n')
      
      cat("### { .mainText}\n\n")
      cat("#### Underexpressed genes \n\n")
      
      cat("Biological pathways affected by underexpressed genes for signature", process, "according to the differential gene expression analysis.\n\n")

      cat('Data obtained from DEG: PATHWAY ANALYSIS UNDEREXPRESSED {.tabset} \n')
      cat('-----------------------------------------------------------------------  \n\n')

      cat("### Gene ontology pathways \n\n")
      cat("```{r table_GO_Under_LPD_", process, ", echo=FALSE}\n", sep = "")
      cat("DT::datatable(params$DEgoUnder_List[[", process, "]], rownames = FALSE, options = list(columnDefs = list(list(className = 'dt-center', targets = '_all'))))\n")
      cat("```\n")
      cat("\n\n")

      cat("### KEGG pathways\n\n")
      cat("```{r table_KEGG_Under_LPD_", process, ", echo=FALSE}\n", sep = "")
      cat("DT::datatable(params$DEkeggUnder_List[[", process, "]], rownames = FALSE, options = list(columnDefs = list(list(className = 'dt-center', targets = '_all'))))\n")
      cat("```\n")
      cat("\n\n")

      cat("### Cancer hallmarks \n\n")
      cat("```{r table_hallmarks_Under_LPD_", process, ", echo=FALSE}\n", sep = "")
      cat("DT::datatable(params$DEHallmarkUnder_List[[", process, "]], rownames = FALSE, options = list(columnDefs = list(list(className = 'dt-center', targets = '_all'))))\n")
      cat("```\n")
      cat("\n\n")

    }

  
  sink()
  
}
