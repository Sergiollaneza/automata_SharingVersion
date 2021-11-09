create.consensus <- function(project, snv_types){

  cat("Starting consensus.... /n")
  consensus <- list()
  for(i in snv_types){
    consensus[[i]] <- read_csv(here("TCGA_results", project, "data", sprintf("snv_%s_data_%s.csv", i, project)))
  }

  indels <- list(consensus[[2]], consensus[[4]])

  # Clean data, check which elements are repeated at least 3 times and filters according to that
  consensus %<>%
    lapply(function(x){x[!names(x) %in% c("X1","ID", "QUAL", "filename")]}) %>%
    lapply(function(x){mutate(x, unique_id = paste(CHROM, POS, REF, ALT, case, sep = "_"))}) %>%
    purrr::reduce(bind_rows)

  count <- table(consensus$unique_id)
  right_ones <- names(count)[as.vector(count) >= 3]
  consensus <- consensus[consensus$unique_id %in% right_ones,]

  # Spreads the info to one row per sample
  filter_info <- consensus %>%
    group_by(unique_id) %>%
    summarize(Info = paste(INFO, collapse = " | "),
              Filter =  paste(FILTER, collapse = " | "))

  consensus %<>% select(-c(FILTER, INFO)) %>%
    distinct %>%
    left_join(filter_info, by = "unique_id")

  # Add the indels
  indels %<>%
    lapply(function(x){filter(x, (nchar(REF) + nchar(ALT)) > 2)}) %>%
    lapply(function(x){x[!names(x) %in% c("X1","ID", "QUAL", "filename")]}) %>%
    lapply(function(x){mutate(x, unique_id = paste(CHROM, POS, REF, ALT, case, sep = "_"))}) %>%
    purrr::reduce(bind_rows)

  count <- table(indels$unique_id)
  right_ones <- names(count)[as.vector(count) >= 2]
  indels <- indels[indels$unique_id %in% right_ones,]

  filter_info <- indels %>%
    group_by(unique_id) %>%
    summarize(Info = paste(INFO, collapse = " | "),
              Filter =  paste(FILTER, collapse = " | "))

  indels %<>% select(-c(FILTER, INFO)) %>%
    distinct %>%
    left_join(filter_info, by = "unique_id")

  # Merge all together and filter by quality
  consensus <- bind_rows(consensus, indels)

  consensus %<>%
    filter(!grepl("normal|germline|homologous|triallelic|Tier1|Tier2|Tier3", Filter)) %>%
    filter(str_extract(Info, "[0-9]{1,3}$") > 40)


  return(consensus)


}















