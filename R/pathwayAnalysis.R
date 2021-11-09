#' Runs the pathway analysis for the given vector names
#'
#' Performs the pathway analysis of the given vector names.
#' @param project .
#' @param automataID .
#' @param process .
#' @param geneNames Genes names to analyse
#' @param contextNames Gene names of the whole dataset (extracted from TCGA)
#' @param positive.l2FC TRUE/FALSE
#' @param type Possible options: "expression" ...
#' @import org.Hs.eg.db
#' @import magrittr
#' @importFrom clusterProfiler bitr enrichKEGG setReadable enricher cnetplot enrichGO dotplot
#' @importFrom msigdbr msigdbr
#' @importFrom dplyr filter select arrange left_join group_by summarise
#' @importFrom tidyr separate_rows
#' @importFrom readr write_csv

runPathwayAnalysis <- function(project, automataID, process, geneNames,
                               contextNames, positive.l2FC, type){

  message("ooo Running pathway analysis for ", process, " ", type)

  if(length(geneNames) > 0){

# Sets the correct name
    if(type == "expression"){
      if(positive.l2FC == TRUE) tag = "overexpressedGenes" else tag = "underexpressedGenes"

    } else if(type == "methylation"){
      if(positive.l2FC == TRUE) tag = "hypermethylatedGenes" else tag = "hypomethylatedGenes"

    } else if(type == "cnv"){
      if(positive.l2FC == TRUE) tag = "amplifiedGenes" else tag = "deletedGenes"

    } else if(type == "snv"){
      tag = "mutatedGenes"

    } else{
      tag = "notDefined"
    }


  # DATA PROCESSING ---------------------------------------------------------


    # Transform gene symbols to entrezID and ensembl which actually work with KEGG and GO
    # For significantly DE genes
    tryCatch(geneNamesTransformed <- bitr(geneNames, fromType = "SYMBOL",
                                 toType=c("ENTREZID", "ENSEMBL"),
                                 OrgDb="org.Hs.eg.db") %>%
      dplyr::filter(is.na(ENTREZID) == FALSE),
      error = function(c){message("WARNING: SYMBOL genes are not mapped...")})

    if(exists("geneNamesTransformed") == FALSE){
      return("Pathway analysis is not available for the selected genes")
    }

    entrezIDselected <- geneNamesTransformed

    # For all of the genes for context
    contextNamesTransformed <- bitr(contextNames, fromType = "SYMBOL",
                                    toType=c("ENTREZID", "ENSEMBL"),
                                    OrgDb="org.Hs.eg.db") %>%
      dplyr::filter(is.na(ENTREZID) == FALSE)

    entrezIDcontext <- contextNamesTransformed


    # The vector can contain ENSEMBL codes and need to be separated
    ensemblNames <- geneNames[substr(geneNames, 1, 4) == "ENSG"]
    ensemblContext <- contextNames[substr(contextNames, 1, 4) == "ENSG"]


    if(length(ensemblNames) > 0){
      if(is.na(ensemblNames[1]) == FALSE){
        tryCatch(ensemblNamesTransformed <- bitr(ensemblNames, fromType = "ENSEMBL",
                                        toType=c("ENTREZID", "SYMBOL"),
                                        OrgDb="org.Hs.eg.db") %>%
          dplyr::filter(is.na(ENTREZID) == FALSE),
          error = function(c){message("WARNING: ENSEMBL genes are not mapped...")})

        # Merges all the genes back
        try(entrezIDselected <- bind_rows(geneNamesTransformed, ensemblNamesTransformed))
      }
    }

    if(length(ensemblContext) > 0){
      if(is.na(ensemblContext[1]) == FALSE){
      ensemblContextTransformed <- bitr(ensemblContext , fromType = "ENSEMBL",
                                      toType=c("ENTREZID", "SYMBOL"),
                                      OrgDb="org.Hs.eg.db") %>%
        dplyr::filter(is.na(ENTREZID) == FALSE)

      # Merges all the genes back
      entrezIDcontext <- bind_rows(contextNamesTransformed, ensemblContextTransformed)

      }
    }


  # PAHYWAY ANALYSIS --------------------------------------------------------


    # Performs KEGG #########
    message("oooo Kegg Analysis for ", process, " ", tag)
    keggAnalysis <- enrichKEGG(gene = entrezIDselected$ENTREZID,
                               organism = "hsa",
                               pvalueCutoff = 0.05,
                               pAdjustMethod = "BH",
                               qvalueCutoff = 0.01)

    if(is.null(keggAnalysis) == FALSE){

      keggAnalysis <- setReadable(keggAnalysis, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

      # Exports KEGG if actually has any result
      if(nrow(keggAnalysis@result) > 0){

        # Exports KEGG table
        write_csv(keggAnalysis@result, here::here("TCGA_results", project, automataID, paste0(type, "_analysis"), "pathway_analysis",
                                                  paste0("keggAnalysis_", tag, "_", process, "_", automataID, ".csv")))


        # Plots only the KEGG if there are significant results
        if(dim(keggAnalysis)[1] > 0){
          # Plots the KEGG analysis
          keggAnalysis.plot <- cnetplot(keggAnalysis, showCategory = 10)


          pdf(here::here("TCGA_results", project, automataID, paste0(type, "_analysis"), "pathway_analysis",
                         paste0("keggAnalysis_", tag, "_", process, "_", automataID, ".pdf")),
              width = 10, height = 10)
          print(keggAnalysis.plot)
          dev.off()

          save(keggAnalysis.plot, file = paste0("TCGA_results/", project, "/", automataID,
                                                "/report_plots/keggAnalysis_Plot_", tag, "_", project, "_", process, "_", automataID, ".Rdata"))
        }
      }

    }


    # Performs GO analysis #########
    message("oooo GO Analysis for ", process, " ", tag)
    ego <- enrichGO(gene = entrezIDselected$ENTREZID,
                    universe = entrezIDcontext$ENTREZID,
                    OrgDb = org.Hs.eg.db,
                    ont = "BP", pAdjustMethod = "BH",
                    qvalueCutoff = 0.01, pvalueCutoff = 0.05,
                    readable = TRUE)

    egoSummary <- as_tibble(ego)


    # Exports GO if actually has any result
    if(nrow(egoSummary) > 0 && !all(is.na(egoSummary))){

      write_csv(egoSummary, here::here("TCGA_results", project, automataID, paste0(type, "_analysis"), "pathway_analysis",
                                       paste0("goAnalysis_", tag, "_", process, "_", automataID, ".csv")))

      # Dotplot
      egoDot.plot <- dotplot(ego, showCategory = 50)

      pdf(here::here("TCGA_results", project, automataID, paste0(type, "_analysis"), "pathway_analysis",
                     paste0("egoDotplot_", tag, "_", process, "_", automataID, ".pdf")),
          width = 20, height = 10)
      print(egoDot.plot)
      dev.off()

      save(egoDot.plot, file = paste0("TCGA_results/", project, "/", automataID,
                                            "/report_plots/egoDot_Plot_", tag, "_", project, "_", process, "_", automataID, ".Rdata"))

      # cnetplot
      cnet.plot <- cnetplot(ego, showCategory = 10, vertex.label.font = 6)

      pdf(here::here("TCGA_results", project, automataID, paste0(type, "_analysis"), "pathway_analysis",
                     paste0("cnetplot_", tag, "_", process, "_", automataID, ".pdf")),
          width = 20, height = 10)
      print(cnet.plot)
      dev.off()

      save(cnet.plot, file = paste0("TCGA_results/", project, "/", automataID,
                                      "/report_plots/cnet_Plot_", tag, "_", project, "_", process, "_", automataID, ".Rdata"))


    }

    # Performs hallmarks analysis #########
    message("oooo Hallmarks Analysis for ", process, " ", tag)
    hallmarksReference <- msigdbr(species = "Homo sapiens", category = "H") %>%
      dplyr::select(gs_name, entrez_gene) %>%
      as_tibble()

    nameReference <- dplyr::select(entrezIDselected, SYMBOL, ENTREZID)

    hallmarks <- enricher(gene = na.omit(entrezIDselected$ENTREZID),
                          universe = entrezIDcontext$ENTREZID,
                          TERM2GENE = hallmarksReference,
                          TERM2NAME = nameReference,
                          pAdjustMethod = "BH", qvalueCutoff = 0.01, pvalueCutoff = 0.05)

    if(is.null(hallmarks) == FALSE){

      hallmarksSummary <- hallmarks@result %>%
        arrange(desc(Count))

      # COnverts the entrezid to hgnc and exports only if there are results
      if(nrow(hallmarksSummary > 0)){

        hallmarksSummary <- hallmarksSummary %>%
          separate_rows(geneID, sep = "/") %>%
          left_join(nameReference, by = c("geneID" = "ENTREZID")) %>%
          group_by(ID, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, Count) %>%
          summarise(entrezID = paste(geneID, collapse = "/"),
                    hgnc_symbol = paste(SYMBOL, collapse = "/")) %>%
          arrange(p.adjust)

        # Exports csv
        write_csv(hallmarksSummary, here::here("TCGA_results", project, automataID, paste0(type, "_analysis"), "pathway_analysis",
                                               paste0("hallmarksAnalysis_", tag, "_", process, "_", automataID, ".csv")))

      }
    }

  }
}
