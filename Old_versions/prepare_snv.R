

prepare.snv <- function(query, snv_files, directory = "GDCdata"){

  source <- ifelse(query$legacy,"legacy","harmonized")
  files <- file.path(query$results[[1]]$project, source,
                     gsub(" ","_",query$results[[1]]$data_category),
                     gsub(" ","_",query$results[[1]]$data_type),
                     gsub(" ","_",query$results[[1]]$file_id),
                     gsub(" ","_",query$results[[1]]$file_name))

  files <- file.path(directory, files)

  # Para todos
  ret <- lapply(files, read.vcfR)
  for(i in 1:length(ret)){
    snps = as.data.frame(ret[[i]]@fix)
    start = str_locate_all(files[i], "/")
    dir = substr(files[i], start = (start[[1]][nrow(start[[1]])]) + 1, stop = nchar(files[1]))
    file_match = snv_files[which(snv_files$file_name == dir),][1,]
    snps$filename = dir
    snps$case = file_match$cases

    if(i == 1){
      snv_data = data.frame(snps)
    } else{
      snv_data = bind_rows(snv_data, snps)
    }
  }
  return(snv_data)
}
