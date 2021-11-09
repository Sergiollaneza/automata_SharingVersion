#' Run LPD stage A
#'
#' Launches the first stage of LPD
#' @param project Name of the project to be run
#' @param automataId Automata ID of the data
#' @param pathway Pathway to the folder that contains the TCGA_results folder. Defaults to current working directory.
#' @importFrom reticulate import source_python
#' @export

runLPD_A <- function(project, automataId, pathway = getwd()){

  message("")
  message("--------------------------------")
  message("~ RUNNING FIRST STAGE OF LPD ~")
  message("--------------------------------")

  Sys.setenv(PATH= paste("/gpfs/software/python/anaconda/4.2/anaconda2/bin",Sys.getenv()["PATH"],sep=";"))
  use_python("/gpfs/software/python/anaconda/4.2/anaconda2/bin/python")

  import("os")
  import("sys")
  import("datetime")
  import("numpy")
  import("csv")




  path <- system.file("python", package = "automata")
  source_python(paste0(path, "/LPD_stageA.py"))

  launch_LPD_A(project, automataId, pathway)

}

