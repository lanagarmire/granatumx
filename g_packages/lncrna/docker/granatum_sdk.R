library(tibble)
library(purrr)
library(dplyr)
library(tidyr)

# .gn and all it's fields are not supposed to be accessed by the user

.gn <- list()
.gn$swd                <- if (Sys.getenv('GRANATUM_SWD') != '') Sys.getenv('GRANATUM_SWD') else '/data'
.gn$args               <- jsonlite::read_json(file.path(.gn$swd, 'args.json'), simplifyVector=T)
.gn$uploaded_files_dir <- file.path(.gn$swd, 'uploaded_files')
.gn$exports_dir        <- file.path(.gn$swd, 'exports')
.gn$exports_anno_file  <- file.path(.gn$swd, 'exports_anno.json')
.gn$imports_dir        <- file.path(.gn$swd, 'imports')
.gn$args_file          <- file.path(.gn$swd, 'args.json')
.gn$results_file       <- file.path(.gn$swd, 'results.json')
.gn$debug_dir          <- file.path(.gn$swd, 'debug')
.gn$dynamic_exports    <- list()
.gn$results            <- list()


# below are all the public methods the user should use

gn_get_uploaded_file_path <- function(inject_into) {
  Sys.glob(file.path(.gn$uploaded_files_dir, inject_into, '*'))[1]
}

gn_get_import <- function(inject_into) {
  jsonlite::read_json(file.path(.gn$imports_dir, inject_into), simplifyVector=T)
}

gn_get_arg <- function(inject_into, default=NULL) {
  arg <- .gn$args[[inject_into]]
  if (is.null(arg)) default else arg
}

# normally you should use this function
gn_export_statically <- function(data, extract_from) {
  jsonlite::write_json(data, file.path(.gn$exports_dir, extract_from), auto_unbox=T)
}

# this is only for those exports that are unkown before running gbox
gn_export_dynamically <- function(data, extract_from, kind=NULL, meta=NULL) {
  .gn$dynamic_exports[[length(.gn$dynamic_exports) + 1]] <<- list(
    extractFrom = extract_from,
    kind = kind,
    meta = meta
  )
  gn_export_statically(data, extract_from)
}

gn_add_ggplot_to_results <- function(g, description=NULL, zoom=2, width=750, height=650, dpi=300) {
  save_filepath <- tempfile()
  ggplot2::ggsave(save_filepath, plot=g, device='png', scale=zoom, width=(width/dpi), height=(height/dpi), dpi=dpi, limitsize=F)
  ggp_raw <- readBin(save_filepath, 'raw', file.info(save_filepath)$size)
  ggp_b64 <- base64enc::base64encode(ggp_raw)
  .gn$results[[length(.gn$results) + 1]] <<- list(
    type = 'png',
    width = width,
    height = height,
    description = description,
    data = ggp_b64
  )
}

gn_add_dataframe_to_results <- function(df, description=NULL, rowname_header=' ') {
  print('rownames(df) =')
  print(rownames(df))
  if (!is.null(rownames(df))) {
    df <- df %>% rownames_to_column(rowname_header)
  }

  table <- df %>%
    as_tibble %>%
    mutate(`__row.numbers__`=row_number()) %>%
    nest(-`__row.numbers__`) %>%
    deframe %>%
    set_names(NULL) %>%
    map(as.list)

  cols <- colnames(df) %>% map(~list(name=.))

  .gn$results[[length(.gn$results) + 1]] <<- list(
    type = 'table',
    description = description,
    cols = cols,
    data = table
  )
}


gn_add_result <- function(data, data_type='json', other_attributes=NULL) {
  .gn$results[[length(.gn$results) + 1]] <<- c(
    list(
      type = data_type,
      data = data
    ),
    other_attributes
  )
}

gn_commit <- function() {
  jsonlite::write_json(.gn$dynamic_exports, .gn$exports_anno_file, auto_unbox=T)
  jsonlite::write_json(.gn$results, .gn$results_file, auto_unbox=T)
}

gn__saveRDS_to_debug <- function(data, filename='debug.RDS') {
  saveRDS(data, file.path(.gn$debug_dir, filename))
}
