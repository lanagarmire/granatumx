library(dplyr)
library(monocle) # Gene filtering with dispersion

source('./granatum_sdk.R') # Uses GranatumX SDK

# Gene filtering function that takes matrix and filters it
# by input parameters - log mean expression and dispersion fit thresholds.
# Along with the result matrix, it should also show a plot showing original
# data among the filtered data (expression vs. dispersion plot).

#
# FOR TESTING
# in_matrix <- in_json$assay$matrix
# in_gids <- in_json$assay$geneIds
# expr_thresh <- -2.3
# disp_thresh <- 1
# temp_plot_filename <- "temp_gene_filtering_plot.png"
#

filter_genes_and_make_plot <- function(in_matrix,
                                       in_gids,
                                       in_sample_ids,
                                       expr_thresh = -2.3,
                                       disp_thresh = 1,
                                       temp_plot_filename = "temp_gene_filtering_plot.png") {
  in_matrix_with_gids <- in_matrix
  rownames(in_matrix_with_gids) <- in_gids
  # The lowerDetectionLimit is lower than default (0.1) or what used before (0.0001)
  cds <- newCellDataSet(in_matrix_with_gids,
                   lowerDetectionLimit = 0,
                   expressionFamily = negbinomial.size())
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds)
  mdt <- dispersionTable(cds)
  mog <- subset(mdt,
           mean_expression >= expr_thresh &
           dispersion_empirical >= disp_thresh * dispersion_fit)$gene_id
  cds <- setOrderingFilter(cds, mog)
  disp_table <- dispersionTable(cds) %>%
    mutate(log_mean_expression      = log(mean_expression)) %>%
    mutate(log_dispersion_empirical = log(dispersion_empirical))

  ordering_genes <-
    row.names(subset(fData(cds), use_for_ordering == TRUE))

  # Return this plot for PNG creation then conversion to JSON
  p <-
    ggplot(disp_table,
           aes(log_mean_expression, log_dispersion_empirical)) +
    geom_point(color = 'darkgrey') +
    geom_line(aes(y = log(dispersion_fit)), color = 'red') +
    geom_point(data = disp_table %>% filter(gene_id %in% ordering_genes),
               color = 'black') +
    labs(x = 'Log Mean Expression', y = 'Log Empirical Dispersion') +
    theme_bw()
  ggsave(filename=temp_plot_filename, plot=p, width=7.5, height=5, dpi=100)
  #p_raw <- readBin(temp_plot_filename, 'raw', file.info(temp_plot_filename)$size)

  fn_output <- list()
  fn_output[["results"]] <- list()
  fn_output[["exports"]] <- list()

  # Results
  # Plot
  # fn_output[["results"]][["gene_filtering_plot"]] <- unbox(base64enc::base64encode(p_raw))
  #fn_output[["results"]][["gene_filtering_plot_raw"]] <- p_raw
  fn_output[["results"]][["gene_filtering_plot_raw"]] <- p

  filter_vec <- rownames(fData(cds))[fData(cds)$use_for_ordering]


  # Genes before and after filtering
  fn_output[["results"]][["numbers_of_genes"]] <-
    #list(list("Before" = unbox(nrow(in_matrix)),
    #     "After" = unbox(length(filter_vec))))
    #list(list("Before" = nrow(in_matrix),
    #     "After" = length(filter_vec)))
    sprintf(paste(c(
    "Numbers of genes:",
    "  - Before = %s",
    "  - After = %s"
    ), collapse="\n"), nrow(in_matrix), length(filter_vec))

  # Filtered geneIDs
  fn_output[["exports"]][["geneIds_filtered"]] <- filter_vec

  # Filtered assay
  tmp_mat <- in_matrix_with_gids[filter_vec, ]
  rownames(tmp_mat) <- NULL

  fn_output[["exports"]][["matrix_filtered"]] <- list(
    matrix = tmp_mat,
    geneIds = filter_vec,
    sampleIds = in_sample_ids
  )

  # Return the assay as a list containing a vector of geneIDs and
  # the reduced matrix, plus before and after filtering numbers
  fn_output
}

# Input parameters/data and calculate DE table

# Import assay and selected option(s)

#
# FOR TESTING
# in_json <- read_json('input_with_genes_de.json', simplifyVector=T)
#

# in_json <- read_json('stdin', simplifyVector=T)

# Explicitly set seed as standard protocol
# rand_seed <- in_json$seed
seed <- gn_get_arg('seed')
set.seed(seed)

assay <- gn_get_import('assay')
expressionThreshold <- gn_get_arg('expressionThreshold')
dispersionThreshold <- gn_get_arg('dispersionThreshold')

# Get exports/results
# output <- filter_genes_and_make_plot(in_json$assay$matrix,
#                                      in_json$assay$geneIds,
#                                      in_json$assay$sampleIds,
#                                      in_json$expressionThreshold,
#                                      in_json$dispersionThreshold)
output <- filter_genes_and_make_plot(assay$matrix,
                                     assay$geneIds,
                                     assay$sampleIds,
                                     expressionThreshold,
                                     dispersionThreshold)

# Exports
#gene_filtered_assay <- list(matrix = output[['exports']][['matrix_filtered']],
#                            geneIds = output[['exports']][['geneIds_filtered']],
#                            sampleIds = assay$sampleIds)
gene_filtered_assay <- assay
gene_filtered_assay$matrix <- output[['exports']][['matrix_filtered']][['matrix']]
gene_filtered_assay$geneIds <- output[['exports']][['matrix_filtered']][['geneIds']]
gene_filtered_assay$sampleIds <- output[['exports']][['matrix_filtered']][['sampleIds']]
gn_export_statically(gene_filtered_assay, 'geneFilteredAssay')

# Results
gn_add_ggplot_to_results(
  output[['results']][['gene_filtering_plot_raw']],
  "Gene filtering plot",
  dpi = 150,
  height = 400
)

gn_add_result(
  output[['results']][['numbers_of_genes']],
  #'text',
  #list(label = 'Number of genes', description = 'Number of genes remaining after filtering.')
  'markdown',
  list(label = 'Number of genes') # the "description" argument is no longer used
)

#
# FOR TESTING
# write_json(output, "temp", simplifyVector=T)
#

# write_json(output, stdout(), simplifyVector=T)
gn_commit()
