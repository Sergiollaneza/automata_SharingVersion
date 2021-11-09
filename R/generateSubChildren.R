generateSubChilds <- function(project, automataID, process, workingPath){
  
  sink(paste0(workingPath, "/TCGA_results/", project, "/", automataID,
              "/Report_RMD/Childs/Signatures/page_signatures_", process, ".Rmd"))
  
  # Creates title
  cat('Signature', process, '{data-navmenu="Signature analysis"}\n')
  cat('=================================================================\n\n')
  
  
  # CREATES CHILD FOR CLINICAL ---------------------------
  if(file.exists(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/clinical_analysis/",
                        "KM_Survival_Overall_Plot_LPD_", process,  "_", automataID, ".pdf"))){
    

    cat('Data obtained from clinical: INFO--CLINICAL\n')
    cat('-----------------------------------------------------------------------\n\n')
    
    cat("### { .mainText}\n\n")
    cat("#### A. Clinical analysis \n\n")
    
    cat("Clinical data was downloaded for each patient but it is only available for tumour samples. This means that if a signature is composed only of normal tissue samples, there is not data available,
        and that the number of samples represented in the tables and plots could not match the actual number of available samples.\n\n")
    
    cat("A Kaplan-Meier estimation was done for all the signatures (Fig. A1) to calculate the survival probability. In addition to that, a Kaplan-Meier comparing the survival probability of the
        samples assigned to the signature", process,  "against all the other samples (Fig. A2). \n\n")

    cat('Data obtained from clinical: KM_PLOTS  \n')
    cat('-----------------------------------------------------------------------\n\n')
    
    cat("### *Figure A1. Kaplan-Meier estimation of survival for all the signatures. *\n\n")
    cat("```{r plot_overall_LPD", process, ", echo=FALSE}\n", sep = "")
    cat("params$KMplot_List[[16]]$plot\n")
    cat("```\n")
    cat("\n\n")
    
    cat("### *Figure A2. Kaplan-Meier estimation of survival for the signature", process, "against all the rest. *\n\n")
    cat("```{r plot_subtype_LPD", process, ", echo=FALSE}\n", sep = "")
    cat("params$KMplot_List[[", process, "]]$plot\n")
    cat("```\n")
    cat("\n\n")
    
    
  }
  
  # CREATES CHILD FOR DIFF EXPRESSED GENES ---------------------------
  if(file.exists(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/report_plots/",
                        "deGenes_Plot_", project, "_", automataID, ".Rdata")) &&
     file.exists(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/expression_analysis/",
                        "DE-Genes_", project, "_LPD_", process, "_", automataID, ".csv"))){


    cat('Data obtained from DEGl: INFO--DIFF EXPRESSED GENES\n')
    cat('-----------------------------------------------------------------------\n\n')

    cat("### { .mainText}\n\n")
    cat("#### B. Differential expressed genes analysis \n\n")

    cat("A differential expression analysis was performed for all the genes between the samples assigned to the signature", process, "against all others (Table B1) using DESeq2. Genes with log fold superior to 1 or
    inferior to -1 were considered significative, with a p.value threshold of 0.05. A volcano plot was made to represent this visually (Figure B1). \n\n")


    cat('Data obtained from DEG: DEG TABLE AND PLOTS  \n')
    cat('-----------------------------------------------------------------------\n\n')

    cat("### *Table B1. Log fold change for each gene representing how differentially expressed they are. Values higher than 1 are considered overexpressed, while values
        lower than -1 are considered underexpressed. *\n\n")
    cat("```{r table_DEG_LPD_", process, ", echo=FALSE}\n", sep = "")
    cat("DT::datatable(params$DEgenes_List[[", process, "]], rownames = FALSE, options = list(columnDefs = list(list(className = 'dt-center', targets = '_all'))))\n")
    cat("```\n")
    cat("\n\n")

    cat("### *Figure B1. Volcano plot for signature", process, " that visually represents the differentially expressed genes. The name of the top 20 differentially expressed genes is shown. *\n\n")
    cat("```{r plot_subtype_Volcano_", process, ", echo=FALSE}\n", sep = "")
    cat("params$volcanoPlot_List[[", process, "]]\n")
    cat("```\n")
    cat("\n\n")


  }

# CREATES CHILD FOR DIFF METHYLATED GENES ---------------------------

  if(file.exists(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/methylation_analysis/",
                        "DM_signifGenes_", project, "_LPD_", process, "_", automataID, ".csv")) &&
     file.exists(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/report_plots/",
                        "DmChromosome_Plot_", project, "_LPD_", process, "_", automataID, ".Rdata")) &&
     file.exists(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/report_plots/",
                        "dmSample_Plot_", project, "_LPD_", process, "_", automataID, ".Rdata"))){


      cat('Data obtained from DMG: INFO--DIFF METHYLATED GENES\n')
      cat('-----------------------------------------------------------------------\n\n')

      cat("### { .mainText}\n\n")
      cat("#### C. Differentially methylated genes \n\n")

      cat("Methylation affects the expression level of the genome. Areas that are hypermethylated have a decrease on their expression level, which is common in cancer for genes with 
          antitumoral activity. The opposite situation is also possible, in which te level of methylation for some genes is decreased, increasing their expression levels. This phenomen is called
          hypomethylation.\n\n")
      
      cat("We have calculated the methylation level for each gene and compared across LPD groups in a differential analysis (Table C1). If a gene had a positive fold change bigger than one was interpreted
          as hypermethylation, and when it was smaller than -1 it was interpreted as hypomethylation. We have also calculated the ratio hypermethylated/hypomethylated genes for each sample (Fig. C1)
          and for each chrosome (Fig.C2). \n\n")

      cat('Data obtained from DMG: DMG TABLE   \n')
      cat('-----------------------------------------------------------------------\n\n')

      cat("### *Table C1. Differentially methylated genes for signature", process, "in comparison to the other groups. meanlogFC bigger than 1 is hypermethylation, 
          while smaller than -1 is hypomethylation *\n\n")
      cat("```{r table_DMG_LPD_", process, ", echo=FALSE}\n", sep = "")
      cat("DT::datatable(params$DMgenes_List[[", process, "]], rownames = FALSE, options = list(columnDefs = list(list(className = 'dt-center', targets = '_all'))))\n")
      cat("```\n")
      cat("\n\n")

      cat('Data obtained from DMG: DMG PLOTS   \n')
      cat('-----------------------------------------------------------------------\n\n')

      cat("### *Figure C1. Ratio hyper/hypomethylation for each sample for signature ", process, ". *\n\n")
      cat("```{r plot_subtype_RatioHyper-hypoSample_", process, ", echo=FALSE}\n", sep = "")
      cat("params$DMsamplePlot_List[[", process, "]]\n")
      cat("```\n")
      cat("\n\n")

      cat("### *Figure C2. Ratio hyper/hypomethylation for each chromosome for signature ", process, ". *\n\n")
      cat("```{r plot_subtype_RatioHyper-hypoChromosome_", process, ", echo=FALSE}\n", sep = "")
      cat("params$DMchromosomePlot_List[[", process, "]]\n")
      cat("```\n")
      cat("\n\n")

  }


# CREATES CHILD FOR DIFF MUTATED GENES ---------------------------

  if(file.exists(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/snv_analysis/",
                        "mutationProportionComparison_", project, "_LPD_", process, "_", automataID, ".csv"))  &&
     file.exists(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/report_plots/",
                        "oncoplot_Plot_", project, "_LPD_", process, "_", automataID, ".Rdata")) &&
     file.exists(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/report_plots/",
                        "snpType_Plot_", project, "_LPD_", process, "_", automataID, ".Rdata")) &&
     file.exists(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/report_plots/",
                        "variantClassification_Plot_", project, "_LPD_", process, "_", automataID, ".Rdata")) &&
     file.exists(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/report_plots/",
                        "variantSample_Plot_", project, "_LPD_", process, "_", automataID, ".Rdata")) &&
     file.exists(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/report_plots/",
                        "variantType_Plot_", project, "_LPD_", process, "_", automataID, ".Rdata"))){


    cat('Data obtained from DMG: INFO--DIFF MUTATED GENES\n')
    cat('-----------------------------------------------------------------------\n\n')

    cat("### { .mainText}\n\n")
    cat("#### D. Single Nucleotide Variants \n\n")

    cat("There are many types of mutations that can happen in the genome and the impact of their effect can vary widely. Mutations can be insertions, deletions of punctual changes of a base 
    for another (Figure D2), and depending on these they can cause missense and nonsese mutations, frame shift or in frame mutations or nonstop mutations (Figure D1). The change of a base is what
    we consider an SNP (Figure D3). We have compared the proportion of mutations present in the signature", process, "against all other to see if there are significant differences (Table D1) and
        in addition to the total number of mutations (Figure D4).. \n\n")

    cat('Data obtained from SNV: SNV TABLE   \n')
    cat('-----------------------------------------------------------------------\n\n')

    cat("### *Table D1. Proportion comparison of the number of mutations for each gene. *\n\n")
    cat("```{r table_mutatedGenes_LPD_", process, ", echo=FALSE}\n", sep = "")
    cat("DT::datatable(params$mutatedGenes_List[[", process, "]], rownames = FALSE, options = list(columnDefs = list(list(className = 'dt-center', targets = '_all'))))\n")
    cat("```\n")
    cat("\n\n")

    cat('Data obtained from DMG: DMG PLOTS   \n')
    cat('-----------------------------------------------------------------------\n\n')

    cat("### *Figure D1. Variant classification", process, ". *\n\n")
    cat("```{r plot_subtype_VariantClassification_", process, ", echo=FALSE}\n", sep = "")
    cat("params$variantClassificationPlot_List[[", process, "]]\n")
    cat("```\n")
    cat("\n\n")

    cat("### *Figure D2. Variant type", process, ". *\n\n")
    cat("```{r plot_subtype_VariantType_", process, ", echo=FALSE}\n", sep = "")
    cat("params$variantType_List[[", process, "]]\n")
    cat("```\n")
    cat("\n\n")

    cat("### *Figure D3. SNP type", process, ". *\n\n")
    cat("```{r plot_subtype_snpType_", process, ", echo=FALSE}\n", sep = "")
    cat("params$snpTypePlot_List[[", process, "]]\n")
    cat("```\n")
    cat("\n\n")

    cat("### *Figure D4. Number variants sample", process, ". *\n\n")
    cat("```{r plot_subtype_numberVariantsSample_", process, ", echo=FALSE}\n", sep = "")
    cat("params$variantSample_List[[", process, "]]\n")
    cat("```\n")
    cat("\n\n")

    # cat('Data obtained from DMG: ONCOPLOT   \n')
    # cat('-----------------------------------------------------------------------\n\n')
    # 
    # cat("### *Figure D6. Oncoplot", process, ". *\n\n")
    # 
    # cat("```{r plot_subtype_Oncoplot_", process, ", echo=FALSE}\n", sep = "")
    # cat(paste0('knitr::include_graphics(paste0(workingPath, "/TCGA_results/", project, "/", automataID,"/snv_analysis/oncoplotContext_Plot_", project, "_LPD_", ', process, ',"_", automataID, ".tiff"))\n' ))
    # cat("```\n")
    # cat("\n\n")
    # 
    # cat("### *Figure D7. Oncoplot", process, ". *\n\n")
    # cat("```{r plot_subtype_Oncoplot2_", process, ", echo=FALSE}\n", sep = "")
    # cat("params$snpTypePlot_List[[", process, "]]\n")
    # cat("```\n")
    # cat("\n\n")

  }

  # CREATES CHILD FOR CNV GENES ---------------------------

  if(file.exists(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/cnv_analysis/",
                        "aberrationsGene_", project, "_LPD_", process, "_", automataID, ".csv")) &&
     file.exists(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/report_plots/",
                        "aberrationSample_plot_", project, "_LPD_", process, "_", automataID, ".Rdata")) &&
     file.exists(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/report_plots/",
                        "aberrationChromosome_Amp_plot_", project, "_LPD_", process, "_", automataID, ".Rdata")) &&
     file.exists(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/report_plots/",
                        "aberrationChromosome_Del_plot_", project, "_LPD_", process, "_", automataID, ".Rdata"))){


    cat('Data obtained from CNV: INFO--DIFF CNV GENES\n')
    cat('-----------------------------------------------------------------------\n\n')

    cat("### { .mainText}\n\n")
    cat("#### E. Copy Number Variation \n\n")

    cat("Copy number variation are cahnges in a region of the DNA that result in gain or loss in copies of such region. It can variate its size from several base pairs to multiple genes. Similarly
        to SNVs, we have compared the proportion of genes affected by amplification (Table E1) and deletions (Table E2), of base pairs (Figure E1), and chromosome (Figure E2 and E3) \n\n")

    cat('Data obtained from CNV: CNV TABLE   \n')
    cat('-----------------------------------------------------------------------\n\n')

    cat("### *Table E1. Proportion comparison of samples affected by amplifications for each gene. *\n\n")
    cat("```{r table_aberratedGenesAmp_LPD_", process, ", echo=FALSE}\n", sep = "")
    cat("DT::datatable(params$aberratedGenesAmp_List[[", process, "]], rownames = FALSE, options = list(columnDefs = list(list(className = 'dt-center', targets = '_all'))))\n")
    cat("```\n")
    cat("\n\n")
    
    cat("### *Table E2. Proportion comparison of samples affected by deletions for each gene. *\n\n")
    cat("```{r table_aberratedGenesDel_LPD_", process, ", echo=FALSE}\n", sep = "")
    cat("DT::datatable(params$aberratedGenesDel_List[[", process, "]], rownames = FALSE, options = list(columnDefs = list(list(className = 'dt-center', targets = '_all'))))\n")
    cat("```\n")
    cat("\n\n")
    
    cat('Data obtained from CNV: CNV PLOTS   \n')
    cat('-----------------------------------------------------------------------\n\n')

    cat("### *Figure E1. Proportion of base pairs affected by amplifications and deletions. *\n\n")
    cat("```{r plot_subtype_aberratedSample_", process, ", echo=FALSE}\n", sep = "")
    cat("params$aberrationSample_List[[", process, "]]\n")
    cat("```\n")
    cat("\n\n")
    
    cat("### *Figure E2. Proportion of bp for each chromosome affected by amplifications. *\n\n")
    cat("```{r plot_subtype_aberratedChromosomeAmp_", process, ", echo=FALSE}\n", sep = "")
    cat("params$aberrationChromosomeAmp_List[[", process, "]]\n")
    cat("```\n")
    cat("\n\n")
    
    cat("### *Figure E3. Proportion of bp for each chromosome affected by deletions. *\n\n")
    cat("```{r plot_subtype_aberratedChromosomeDel_", process, ", echo=FALSE}\n", sep = "")
    cat("params$aberrationChromosomeDel_List[[", process, "]]\n")
    cat("```\n")
    cat("\n\n")


  }

  # CREATES CHILD FOR INTERSECTION---------------------------

  if(file.exists(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/intersect/",
                        "overexpressedGenesSummary", "_LPD_", process, "_", automataID, ".csv")) &&
     file.exists(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/intersect/",
                        "underexpressedGenesSummary", "_LPD_", process, "_", automataID, ".csv"))){


    cat('Data obtained from INTERSECT: INFO--intersect GENES\n')
    cat('-----------------------------------------------------------------------\n\n')

    cat("### { .mainText}\n\n")
    cat("#### F. Intersect \n\n")

    cat("An intersect of the genes that presented different type of alterations was performed for overexpressed genes (Table F1 and Figure F1),
        and for underexpressed genes (Table F2 and Figure F2). It was considered that genes with mutations can only cause underexpression.\n\n")


    cat('Data obtained from Intersect: Intersect TABLE   \n')
    cat('-----------------------------------------------------------------------\n\n')

    cat("### *Table F1. Intersect of overexpressed genes. *\n\n")
    cat("```{r table_interOverGenes_LPD_", process, ", echo=FALSE}\n", sep = "")
    cat("DT::datatable(params$intersectOverTable_List[[", process, "]], rownames = FALSE, options = list(columnDefs = list(list(className = 'dt-center', targets = '_all'))))\n")
    cat("```\n")
    cat("\n\n")

    cat("### *Figure F1. Venn diagram of overexpressed genes. *\n\n")
    cat("```{r plot_subtype_aberratedSample2_", process, ", echo=FALSE}\n", sep = "")
    cat('plot(params$vennOver_List[[', process, ']], doWeights = FALSE)\n')
    cat("```\n")
    cat("\n\n")

    cat('Data obtained from Intersect: Intersect TABLE2   \n')
    cat('-----------------------------------------------------------------------\n\n')

    cat("### *Table F2. Intersect of underexpressed genes. *\n\n")
    cat("```{r table_interUnderGenes_LPD_", process, ", echo=FALSE}\n", sep = "")
    cat("DT::datatable(params$intersectUnderTable_List[[", process, "]], rownames = FALSE, options = list(columnDefs = list(list(className = 'dt-center', targets = '_all'))))\n")
    cat("```\n")
    cat("\n\n")

    cat("### *Figure F2. Venn diagram of overexpressed genes. *\n\n")
    cat("```{r plot_subtype_aberratedSample3_", process, ", echo=FALSE}\n", sep = "")
    cat('plot(params$vennUnder_List[[', process, ']], type = "ellipses")\n')
    cat("```\n")
    cat("\n\n")


  }

  
  
  sink()
  
}
