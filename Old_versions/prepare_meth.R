prepare.meth <- function(query, meth_files, directory = "GDCdata"){

  source <- ifelse(query$legacy,"legacy","harmonized")
  files <- file.path(query$results[[1]]$project, source,
                     gsub(" ","_",query$results[[1]]$data_category),
                     gsub(" ","_",query$results[[1]]$data_type),
                     gsub(" ","_",query$results[[1]]$file_id),
                     gsub(" ","_",query$results[[1]]$file_name))



  files <- file.path(directory, files)

  for(i in 1:length(files)){
    meth <- read_delim(files[[i]], delim = "\t")
    start = str_locate_all(files[i], "/")
    dir = substr(files[i], start = (start[[1]][nrow(start[[1]])]) + 1, stop = nchar(files[1]))
    file_match = meth_files[which(meth_files$file_name == dir),][1,]
    meth$filename <- dir
    meth$case <- file_match$cases

    if(i == 1){
      meth_data = data.frame(meth)
    } else{
      meth_data = bind_rows(meth_data, meth)
    }
  }
  return(meth_data)
}
