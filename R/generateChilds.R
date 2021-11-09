createChildsRMD <- function(project, automataID, workingPath, LPDvector){


  # DATA EXTRACTED FROM BARCODES --------------------------------------------
  sink(paste0(workingPath, "/TCGA_results/", project, "/", automataID,
              "/Report_RMD/Childs/Barcodes/page_barcodes.Rmd"))

  cat('\n```{r functions_page_barcodes, include=FALSE}\n')
  cat('n_primary <- nrow(dplyr::filter(params$barcodes, Type == "Primary Solid Tumor"))\n')
  cat('n_metastatic <- nrow(dplyr::filter(params$barcodes, Type == "Metastatic"))\n')
  cat('n_normal <- nrow(dplyr::filter(params$barcodes, Type == "Solid Tissue Normal"))\n')
  cat('```\n\n\n')

  cat("Samples\n")
  cat("===========================================\n\n")

  cat('Data obtained from Barcodes: INFO\n')
  cat('-----------------------------------------------------------------------\n\n')

  cat("### { .mainText}\n\n")
  cat("#### Samples info\n\n")
  cat("A total of **`r n_samples` samples for `r params$project`** with available **`r params$method`**
      data had been found (Table 1). Of the `r n_samples` samples, **`r n_patients` belong to
      different patients**.\n")
  cat("According to the sample type, there are a total of `r n_primary` primary tumour samples (`r round(100*(n_primary/n_samples), 2)`% of the total),
      `r n_metastatic` metastatic samples (`r round(100*(n_metastatic/n_samples), 2)`% of the total)
      and `r n_normal` normal tissue samples (`r round(100*(n_normal/n_samples), 2)`% of the total).\n\n")

  cat('Data obtained from Barcodes: TABLE BARCODES\n')
  cat('-----------------------------------------------------------------------\n\n')

  cat("### *Table 1: Participant and sample barcodes for `r params$project` with available `r params$method` data in addition to the sample type and the assigned Automata ID.*\n\n")

  cat("```{r table_barcodes, echo=FALSE}\n")
  cat("DT::datatable(params$barcodes[,-1], options = list(columnDefs = list(list(className = 'dt-center', targets = '_all'))))\n")
  cat("```\n")
  cat("\n")

  sink()


  # DATA EXTRACTED FROM TOP GENES --------------------------------------------

  sink(paste0(workingPath, "/TCGA_results/", project, "/", automataID,
              "/Report_RMD/Childs/Genes/page_genes.Rmd"))

  cat("Genes\n")
  cat("===========================================\n\n")

  cat('Data obtained from Genes: INFO\n')
  cat('-----------------------------------------------------------------------\n\n')

  cat("### { .mainText}\n\n")
  cat("#### Genes info \n\n")
  cat("**Expression `r params$method` was downloaded for each of the `r n_patients` cases**. It may happen that the same case has multiple samples from different
      vials or plates, and thus, **duplicated samples are removed**. If a case has multiple vials, the ones labelled as 'A' have priority over 'B', 'C'
      and subsequent. The reason behind this is that vials B and C are FFPE, meaning that their results are less reliable. On the other hand, if a case has
      multiple plates, the one with the highest lexicographical value is the one chosen as it is assumed that the sample had to be run again due to a mistake.\n\n")

  cat("**Expression data needs to be normalised to run LPD on it**. As in this case is `r params$method`, the normalisation process consist on ")
  cat('`r if(params$method == "RNA-seq"){"using the package DESeq2 which runs a variance stabilizing transformation process in combination to log2."}`')
  cat('`r if(params$method == "microarray"){"taking the minimum value of each column (gene) and sum it up to all the values of the column (samples), so negative values
      become positive and the distance between the values is still the same."}`\n\n')

  cat("Once the expression data is normalised, **the 500 genes with the highest expression are selected to run LPD on them** (Table 1).
      This is due to two reasons, (i) LPD is a complex process and the time it takes to run is dramatically increased for each gene and sample, therefore
      running the whole human genome could possibly take months if not years, and (ii) previous literature has shown that the top 500 genes
      are representative enough of the variance across the pool of samples.\n\n")

  cat('Data obtained from Genes: TABLE TOP 500\n')
  cat('-----------------------------------------------------------------------\n\n')

  cat("### *Table 1: Top 500 genes with the highest variation across all the samples for `r params$project`.*\n\n")
  cat("```{r table_top500genes, echo=FALSE}\n")
  cat("DT::datatable(params$genes, rownames = FALSE, colnames = rep('', ncol(params$genes)), options=list(ordering=F))\n")
  cat("```\n")
  cat("\n")




  sink()


  # DATA EXTRACTED FROM LPD PARAMETERS  --------------------------------------------

  sink(paste0(workingPath, "/TCGA_results/", project, "/", automataID,
              "/Report_RMD/Childs/LPD_parameters/page_LPD.Rmd"))

  cat("LPD parameters\n")
  cat("===========================================\n\n")

  cat('Data obtained from LPD parameters: INFO1--ESTIMATNG PARAMETERS\n')
  cat('-----------------------------------------------------------------------\n\n')

  cat("### { .mainText}\n\n")
  cat("#### Latent process decomposition \n\n")

  cat("**Latent process decomposition (LPD)** is a hierarchical Bayesian model for unsupervised clustering developed by
      [*Rogers et al.*](https://www.researchgate.net/publication/6751775_The_Latent_Process_Decomposition_of_cDNA_Microarray_Data_Sets).
      This algorithm is able to **analyse the gene expression levels from a pool of samples and classify them into several groups** (also known as clusters).
      However, LPD is a fuzzy clustering technique, which means that **one sample can belong to different groups at the same time but with
      different membership levels**.\n\n Imagine that we have a group of samples in which we have detected two groups; a traditional approach would be to classify
      each sample into either group 1 or group 2, but what LPD does is to estimate how much each sample is represented by each group, so it is possible to have a sample
      that is 60% group 1 and 40% group 2, or 100% group 1 - 0% group 2, or any possible combination that sums up 100%. The use of a fuzzy clustering approach
      seems to do a better job in representing the high genomic heterogeneity present in cancer cells, as previous research has proved for
      [prostate cancer](https://pubmed.ncbi.nlm.nih.gov/28753852/).\n\n")

  cat("Since the words *group* and *cluster* have the connotation of a close set of elements, and that does not match with what LPD does, we
      use the term *process*. We also use the term *gamma values* as a unit to measure the level of membership of each sample to each process.\n\n")

  cat("-----------------------------\n\n")

  # cat("&nbsp;\n\n")

  cat("#### Estimating parameters \n\n")

  cat("LPD is an unsupervised machine learning algorithm, therefore it does not require any type of extra data besides the expression levels of the 500 genes with
      the highest variance. Due to this, **LPD has to estimate several parameters to run the analysis from which the most important ones are the number of processes and the
      sigma value**. \n\n")

  cat("In order to do this, **LPD is run multiple times using different combinations of these parameters and the best three are picked** (Table 1). The range of number
      of processes goes from 2 to 15, while the sigma value goes from -0.001 to -1.5 in regular intervals. Every combination of these values is repeated 10 times
      for cross-validation and their fitness is measured by calculating the log-likelihood of the combination mean across all repetitions (Figure 1).\n\n")

  cat("For clinical purposes, it is better to have a low number of subtypes that are very different to each other. We achieve this by, instead of choosing
      the combination with the lowest log-likelihood, choosing a combination within the same range of standard error and with the lowest possible number of
      subtypes. Then, we choose the combinations with &pm;1 process for later validation.\n\n")


  cat('Data obtained from LPD: LIKELIHOOD + BEST COMBOS  \n')
  cat('-----------------------------------------------------------------------\n\n')

  cat("### *Table 1. Three best combinations of number of processes and sigma value.*\n\n")
  cat("```{r table_processAndSigma, echo=FALSE}\n")
  cat("DT::datatable(params$processAndSigmas, rownames = FALSE, options = list(columnDefs = list(list(className = 'dt-center', targets = '_all')), dom = 't'))\n")
  cat("```\n")
  cat("\n")

  cat("### *Figure 1. Mean log-likelihood for each combination. *\n\n")
  cat("```{r plot_likelihood, echo=FALSE}\n")
  cat("params$likelihood_plot\n")
  cat("```\n")
  cat("\n\n")

  cat('Data obtained from LPD parameters: INFO2--CORRELATION\n')
  cat('-----------------------------------------------------------------------\n\n')

  cat("### { .mainText}\n\n")
  cat("#### Selecting the optimal run \n\n")

  cat("**LPD is run again 100 times for each combination** to have more accurate results. Then, the one hundred runs per combinations are compared against each
      other, and the one closest to their mean is selected as the representative of the combination. As stated before, it is better to have subtypes that are
      very differentiated between each other, so we compare the internal correlation levels of each run and the ones with the lowest is picked as the optimal
      run of the whole dataset (Table 2).\n\n")

  cat("For `r params$project` using `r params$method`, **the best-found combination was `r params$bestCombination$n_process` processes with a sigma value
      of `r params$bestCombination$sigma`**. The distribution of the gamma values of this combination can be seen in Figure 2, this depicts how much each
      process is representing the genetic variation of each sample. According to in which process each sample is more represented, they are assigned to different
      groups.\n\n")


    cat('Data obtained from CORRELATION TABLE  \n')
  cat('-----------------------------------------------------------------------\n\n')

  cat("### *Table 2. Internal correlations of each combination within its processes. The correlation was calculated using Pearson correlation coefficient.*\n\n")
  cat("```{r table_correlations, echo=FALSE}\n")
  cat("DT::datatable(params$correlations, rownames = FALSE, options = list(columnDefs = list(list(className = 'dt-center', targets = '_all'))))\n")
  cat("```\n")
  cat("\n")

  cat("### *Figure 2. Gamma values of each sample for each process, the colour indicates the different groups to which samples are assigned to. *\n\n")
  cat("```{r plot_likelihood2, echo=FALSE}\n")
  cat("params$gammaSample\n")
  cat("```\n")
  cat("\n\n")

  sink()

  
  # DATA EXTRACTED FROM SIGNATURES --------------------------------------------
  
  sink(paste0(workingPath, "/TCGA_results/", project, "/", automataID,
              "/Report_RMD/Childs/Signatures/page_signatures.Rmd"))
  
  # Creates the frame for the subtypes child only if they exist
  for(process in 1:length(LPDvector)){
    
    if(LPDvector[process] == TRUE){
      cat('```{r signatures_subchild, child="', workingPath, "/TCGA_results/", project, "/", automataID, "/Report_RMD/", 'Childs/Signatures/page_signatures_', process, '.Rmd"}\n', sep = "")
      cat('```\n')
      cat('\n')
    }
    
  }

  sink()
  
  
  # DATA EXTRACTED FROM OVERALL --------------------------------------------
  sink(paste0(workingPath, "/TCGA_results/", project, "/", automataID,
              "/Report_RMD/Childs/Overall/page_overall.Rmd"))
  
  cat("Cancer analysis \n")
  cat("===========================================\n\n")
  
  
  if(file.exists(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/batch_effect/",
                        "sampleTypeProportions_", project, "_", automataID, ".csv"))){
    
    
    cat('Description of BATCH EFFECT\n')
    cat('-----------------------------------------------------------------------\n\n')
    
    cat("### { .mainText}\n\n")
    cat("#### A. Batch effect\n\n")
    cat("Batch effect is a phenomem that occurs when non-biological factors affect the outcome of an analysis, leading to erroneous conclusions. One of the
        posible reasons for this to happen is that the samples from the TCGA have a mix of tumour and not tumour samples which creates very different
        expression profiles. LPD can catch that as a possible cancer subtype when in reality is just the differences between the tissues. Due to this,
        the number of samples for each type was calculated and a Chi Square test was used to see if there were significant differences in the distribution
        of each sample type (Table A1).\n\n")

    cat('Visual representation of BATCH EFFECT\n')
    cat('-----------------------------------------------------------------------\n\n')
    
    cat("### *Table A1: Quantity of each sample type for each subtype to detect possible batch effect. 01 = Tumour sample, 02 = Recurrent Solid Tumour. A = Frozen tissue, B and C are FFPE.*\n\n")
    
    cat("```{r table_batchEffect, echo=FALSE}\n")
    cat("DT::datatable(params$batchEffect, options = list(columnDefs = list(list(className = 'dt-center', targets = '_all'))))\n")
    cat("```\n")
    cat("\n")
    
    
  }
  
  if(file.exists(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/clinical_analysis/",
                        "raceProportions_", project, "_", automataID, ".csv"))){
    
    cat('Description of RACE PROPORTION\n')
    cat('-----------------------------------------------------------------------\n\n')
    
    cat("### { .mainText}\n\n")
    cat("#### B. Race proportion\n\n")
    cat("Races have differences in the gene expression that can cause LPD to group together all the members
        from a specific race. To check if that is happening, a Chi Square Proportion analysis
        was performed to compare across all signatures if they have the same proportion of race 
        representation (Table B1).\n\n")
    
    cat('Data obtained from Overall: TABLE RACE PROPORTION\n')
    cat('-----------------------------------------------------------------------\n\n')
    
    cat("### *Table B1: Comparison of race proportion for each signature.*\n\n")
    
    cat("```{r table_raceProportion, echo=FALSE}\n")
    cat("DT::datatable(params$raceProportion, options = list(columnDefs = list(list(className = 'dt-center', targets = '_all'))))\n")
    cat("```\n")
    cat("\n")
    
    
  }
  
  
  sink()
  
  
  # DATA EXTRACTED FROM PATHWAYS --------------------------------------------
  
  sink(paste0(workingPath, "/TCGA_results/", project, "/", automataID,
              "/Report_RMD/Childs/Signatures/page_pathways.Rmd"))
  
  # Creates the frame for the subtypes child only if they exist
  for(process in 1:length(LPDvector)){
    
    if(LPDvector[process] == TRUE){
      cat('```{r signatures_subpathway, child="', workingPath, "/TCGA_results/", project, "/", automataID, "/Report_RMD/", 'Childs/Signatures/page_pathways_', process, '.Rmd"}\n', sep = "")
      cat('```\n')
      cat('\n')
    }
    
  }
  
  sink()
  

}
