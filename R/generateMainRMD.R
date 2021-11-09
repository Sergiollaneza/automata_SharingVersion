createMainRMD <- function(project, automataID, workingPath){

  # Strings that need to go quoted and cause double quote
  dayFormat <- format(Sys.Date(), format = "%d/%m/%Y")

  # Creates the script for the main document
  sink(paste0(workingPath, "/TCGA_results/", project, "/", automataID, "/Report_RMD/mainDocument.Rmd"))


# YAML field --------------------------------------------------------------


  cat('---\n')
  cat('title: "Automata report 3.0"\n')
  cat('author: "Sergio Llaneza-Lago"\n')
  cat('date: ', dayFormat, '\n', sep = "")
  cat('description: "Main page for Automata report"\n')
  cat('output: \n')
  cat(' flexdashboard::flex_dashboard:\n')
  cat('  vertical_layout: scroll\n')
  cat('  theme: flatly\n')
  cat('  orientation: rows\n')
  cat('  source_code: https://github.com/UEA-Cancer-Genetics-Lab/Automata\n')
  cat('params:\n')
  cat(' project: NULL\n')
  cat(' automataID: NULL\n')
  cat(' method: NULL\n')
  cat(' path: NULL\n')
  cat(' show_barcodes: FALSE\n')
  cat(' barcodes: NULL\n')
  cat(' show_genes: FALSE\n')
  cat(' genes: NULL\n')
  cat(' show_LPD: FALSE\n')
  cat(' processAndSigmas: NULL\n')
  cat(' likelihood_plot: NULL\n')
  cat(' bestCombination: NULL\n')
  cat(' correlations: NULL\n')
  cat(' sampleAssigned: NULL\n')
  cat(' gammaSample_plot: NULL\n')
  cat(' KMplot_List: NULL\n')
  cat(' DEgenes_List: NULL\n')
  cat(' volcanoPlot_List: NULL\n')
  cat(' DMgenes_List: NULL\n')
  cat(' DMsamplePlot_List: NULL\n')
  cat(' DMchromosomePlot_List: NULL\n')
  cat(' mutatedGenes_List: NULL\n')
  cat(' oncoplot_List: NULL\n')
  cat(' snpTypePlot_List: NULL\n')
  cat(' variantClassificationPlot_List: NULL\n')
  cat(' variantSample_List: NULL\n')
  cat(' variantType_List: NULL\n')
  cat(' aberratedGenesAmp_List: NULL\n')
  cat(' aberratedGenesDel_List: NULL\n')
  cat(' aberrationSample_List: NULL\n')
  cat(' aberrationChromosomeAmp_List: NULL\n')
  cat(' aberrationChromosomeDel_List: NULL\n')
  cat(' intersectOverTable_List: NULL\n')
  cat(' intersectUnderTable_List: NULL\n')
  cat(' batchEffect: NULL\n')
  cat(' raceProportion: NULL\n')
  cat(' DEgoOver_List: NULL\n')
  cat(' DEgoUnder_List: NULL\n')
  cat(' DEkeggOver_List: NULL\n')
  cat(' DEkeggUnder_List: NULL\n')
  cat(' DEHallmarkOver_List: NULL\n')
  cat(' DEHallmarkUnder_List: NULL\n')
  cat(' vennUnder_List: NULL\n')
  cat(' vennOver_List: NULL\n')
  cat('\n')
  cat('---\n')
  cat('\n')

  # cat('\n')
  # cat('<style type="text/css">\n\n')
  # cat('.chart-stage {\n')
  # cat('  font-size: 18px;\n')
  # cat("}\n\n")

  cat('\n')
  cat('<style type="text/css">\n\n')
  # Changes the text that appears explaining everything
  cat('.mainText {\n')
  cat('  font-size: 17px;\n')
  cat("}\n\n")

  # Creates a fake header -- EXPERIMENTAL
  cat('.blackFakeHeader {\n')
  cat('  font-size: 22px;\n')
  cat('  color: black;\n')
  cat('  float: left;\n')
  cat("}\n\n")

  # Changes the size of the description of tables/figures
  cat('.chart-title {\n')
  cat('  font-size: 16px;\n')
  cat("}\n\n")

  # Changes the colour and size of the header level 4
  cat('h4 {\n')
  cat('  font-size: 22px;\n')
  cat('  color: #00A989;\n')
  cat("}\n\n")

  cat('</style>\n\n')


# OVERVIEW PAGE -----------------------------------------------------------


  cat('<!-- ######### 1. OVERVIEW PAGE ######### -->\n\n')

  cat('Overview\n')
  cat('=========================================================================================\n\n')

  # Barcodes overview
  cat('<!-- 1. 1. Barcodes overview child-->\n\n')
  cat('```{r barcodes_overview_child, child="Childs/Barcodes/overview_barcodes.Rmd", eval = params$show_barcodes}\n')
  cat('```\n')
  cat('\n')
  # LPD overview
  cat('<!-- 1. 2. LPD overview child-->\n\n')
  cat('```{r lpd_overview_child, child="Childs/LPD_parameters/overview_LPD.Rmd", eval = params$show_LPD}\n')
  cat('```\n')
  cat('\n')


# BARCODES PAGE -----------------------------------------------------------

  cat('<!-- 2. Barcodes page -->\n\n')
  cat('```{r barcodes_page_child, child="Childs/Barcodes/page_barcodes.Rmd", eval = params$show_barcodes}\n')
  cat('```\n')
  cat('\n')



# TOP GENES PAGE -----------------------------------------------------------

  cat('<!-- 3. Top genes page -->\n\n')
  cat('```{r genes_page_child, child="Childs/Genes/page_genes.Rmd", eval = params$show_genes}\n')
  cat('```\n')
  cat('\n')

# LPD PARAMETERS PAGE -----------------------------------------------------------

  cat('<!-- 4. LPD parameters page -->\n\n')
  cat('```{r lpd_page_child, child="Childs/LPD_parameters/page_LPD.Rmd", eval = params$show_LPD}\n')
  cat('```\n')
  cat('\n')

# CANCER ANALYSIS PAGE -----------------------------------------------------------

  cat('<!-- 5. Cancer analysis page -->\n\n')
  cat('```{r cancer_page_child, child="Childs/Overall/page_overall.Rmd", eval = params$show_LPD}\n')
  cat('```\n')
  cat('\n')

# SIGNATURE ANALYSIS PAGES -----------------------------------------------------------

  cat('<!-- 6. Signatury analysis page -->\n\n')
  cat('```{r signatures_page_child, child="Childs/Signatures/page_signatures.Rmd", eval = params$show_LPD}\n')
  cat('```\n')
  cat('\n')

  # SIGNATURE ANALYSIS PAGES -----------------------------------------------------------

  cat('<!-- 6. Pathway analysis page -->\n\n')
  cat('```{r signatures_page_child, child="Childs/Signatures/page_pathways.Rmd", eval = params$show_LPD}\n')
  cat('```\n')
  cat('\n')



  sink()

}
