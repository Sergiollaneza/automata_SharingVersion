#' Function to remove the duplicated plates
#'
#' Removes duplicated samples with different plates or/and vials.
#' @param project Name of a TCGA project (eg. TCGA-PRAD or TCGA-BRCA).
#' @param automataID ID assigned to the analysis so the files from different analysis are not mixed.
#' @param microarray Wether the expression data comes from RNA-seq counts or microarray levels
#' @return A matrix without duplicates
#' @importFrom tibble column_to_rownames rownames_to_column
#' @importFrom tidyr separate
#' @importFrom dplyr group_by arrange slice ungroup select
#' @import magrittr
#' @export

cleanMatrix <- function(project, automataID, microarray = FALSE){

  type <- ifelse(microarray == TRUE, "microarray", "RNA-seq")

  # Reads the matrix
  expressionMatrix <- read.csv(here::here("TCGA_results", project, automataID, "TCGA_data",
                                          sprintf("Expression_matrix_%s_%s_%s.csv", type, project, automataID)))

  if(type == "microarray"){

    expressionMatrix <- column_to_rownames(expressionMatrix, "X")

  } else {

    expressionMatrix <- column_to_rownames(expressionMatrix, "X1")
  }


  # Transposes and separates hte barcode in different columns
  expressionMatrix_t <- as.data.frame(t(expressionMatrix)) %>%
    rownames_to_column("Barcode") %>%
    separate(Barcode, c("Project", "TSS", "Participant", "Sample_Vial", "Portion_Analyte", "Plate", "Center"), "\\.", remove = FALSE) %>%
    separate(Sample_Vial, c("Sample", "Vial"), "(?<=[0-9])(?=[A-Z])") %>%
    separate(Portion_Analyte, c("Portion", "Analyte"), "(?<=[0-9])(?=[A-Z])")

  # Removes rows with same participant and sample but different vial or plate
  # Gives preference to A vials as usually B and C are FFPE
  # Gives preference to the plates with the highest lexicographical sort value
  expressionMatrix_t_clean <- expressionMatrix_t %>%
    group_by(TSS, Participant, Sample) %>%
    arrange(Vial, desc(Plate)) %>%
    dplyr::slice(1)

  # Transposes again and saves it
  expressionMatrix_clean <- expressionMatrix_t_clean %>%
    ungroup() %>%
    select(-c(Project, TSS, Participant, Sample, Vial, Portion, Analyte, Plate, Center)) %>%
    t()

  colnames(expressionMatrix_clean) <- expressionMatrix_clean[1,]

  write.csv(expressionMatrix_clean[2:nrow(expressionMatrix_clean),], here::here("TCGA_results", project, automataID, "processed_expression",
                                                                                sprintf("clean_matrix_%s_%s_%s.csv", type, project, automataID)))

}


#' Normalise expression matrix
#'
#' Applies DESeq2 to an RNA-seq / microarray expression matrix.
#' @param project Name of a TCGA project (eg. TCGA-PRAD or TCGA-BRCA).
#' @param automataID ID assigned to the analysis so the files from different analysis are not mixed.
#' @param microarray Wether the expression data comes from RNA-seq counts or microarray levels
#' @return   A dataframe with samples in the colums and gene codes in the rows containing the
#' normalised expression count values
#' @examples
#' normaliseExpression(expressionMatrix)
#' @importFrom DESeq2 DESeqDataSetFromMatrix vst
#' @export

normaliseExpression <- function(project, automataID, microarray = FALSE){

  message("")
  message("----------------------------------")
  message("~ NORMALISING EXPRESSION MATRIX ~")
  message("----------------------------------")

  # Creates folders to save the data
  message("o Creating folders...")
  createTCGAfolder(project, automataID)

  if(microarray == TRUE){

    message("o Importing microarray data")

    type = "microarray"
    expressionMatrix <- read.csv(here::here("TCGA_results", project, automataID, "processed_expression",
                                            sprintf("clean_matrix_%s_%s_%s.csv", type, project, automataID)), stringsAsFactors = FALSE) %>%
      distinct(X, .keep_all = TRUE)

    rownames(expressionMatrix) <- expressionMatrix$X
    expressionMatrix <- as.matrix(expressionMatrix[,-1])

    message("o Adding the minimum value to all the values to remove negatives")

    # Removes negative values
    normal_matrix <- apply(expressionMatrix, 2, function(x){
      minimum <- min(as.numeric(x), na.rm = TRUE)
      sum <- as.numeric(x) + abs(as.numeric(minimum))
      sum
    })

    rownames(normal_matrix) <- rownames(expressionMatrix)
    normal_matrix <- na.omit(normal_matrix)

  } else{

    message("o Importing RNA-seq data")

    type = "RNA-seq"
    expressionMatrix <- read.csv(here::here("TCGA_results", project, automataID, "processed_expression",
                                            sprintf("clean_matrix_%s_%s_%s.csv", type, project, automataID)), row.names = 1)


    col_data <- data.frame(cbind(colnames(expressionMatrix), condition = "cancer"))
    message("o Applying DESeq")
    des <- DESeqDataSetFromMatrix(countData = na.omit(expressionMatrix), colData = col_data, design = ~ 1)
    vst <- vst(des)
    message("o Generating normalised data frame")
    normal_matrix <- assay(vst)

  }

  message("o Exporting results")

  write.csv(normal_matrix, here::here("TCGA_results", project, automataID, "processed_expression",
                                      sprintf("normal_matrix_%s_%s_%s.csv", type, project, automataID)))
  message(paste("Normal matrix has been saved in ",
                here::here("TCGA_results", project, automataID, "processed_expression",
                           sprintf("normal_matrix_%s_%s_%s.csv", type, project, automataID))))
}


#' Select the genes with the highest variance
#'
#' Order the expression matrix according to the variance of each column and takes the minimum
#' number of needed columns to reach the given number of unique genes.
#' @param project Name of a TCGA project (eg. TCGA-PRAD or TCGA-BRCA).
#' @param automataID ID assigned to the analysis so the files from different analysis are not mixed.
#' @param uniqueGenes Number of desired unique genes in the output.
#' @return Matrix containing the columns with the highest variance. The number of columns depends on the
#' number of uniqueGenes.
#' @export

topGenes <- function(project, automataID, microarray = FALSE, uniqueGenes = 500){

  message("")
  message("-------------------------")
  message("~ SELECTING TOP GENES ~")
  message("-------------------------")

  type <- ifelse(microarray == TRUE,  "microarray", "RNA-seq")

  expressionMatrix <- read.csv(here::here("TCGA_results", project, automataID, "processed_expression",
                                          sprintf("normal_matrix_%s_%s_%s.csv", type, project, automataID)), row.names = 1)

  expressionMatrix <- t(expressionMatrix)

  # Calculate the variance of each column (columns are genes
  message("o Calculating variance of each column...")
  columnVariance <- apply(expressionMatrix, 2, var)
  orderedGenes <- expressionMatrix[, order(columnVariance , decreasing = TRUE)]

  # Take the minimum number of required genes to reach the desired quantitiy
  message("o Selecting top genes...")
  genes <- colnames(orderedGenes)

  i <- uniqueGenes
  while(length(unique(genes[1:i])) != uniqueGenes){
    i <- i + 1
  }

  topGenes <- orderedGenes[,1:uniqueGenes]

  write.csv(topGenes, here::here("TCGA_results", project, automataID, "processed_expression",
                                 sprintf("topGenes_%s_%s.csv", project, automataID)))
  message(paste("Top genes had been saved in ",
                here::here("TCGA_results", project, automataID, "processed_expression",
                           sprintf("topGenes_%s_%s.csv", project, automataID))))

}
