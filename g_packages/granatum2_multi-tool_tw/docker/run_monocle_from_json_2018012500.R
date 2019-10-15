# twolf 20180125 from myscript_2018012504.R
suppressPackageStartupMessages({
  library(monocle)
  library(jsonlite)
  library(ggplot2)
  library(base64enc)
})
sink("/dev/null")
in_json <- read_json('stdin', simplifyVector=T)
expr_matrix = in_json$assay
#run_monocle <- function(expr_matrix) {
  # Expects raw counts
  monocle_sample_sheet <- data.frame(cell_number=c(1:dim(expr_matrix)[2]))
  rownames(monocle_sample_sheet) <- c(1:dim(expr_matrix)[2])
  # The namining "gene_short_name" is necessary
  monocle_gene_annotation <- data.frame(gene_short_name=c(1:dim(expr_matrix)[1]))
  rownames(monocle_gene_annotation) <- c(1:dim(expr_matrix)[1])
  colnames(expr_matrix) <- monocle_sample_sheet$cell_number
  rownames(expr_matrix) <- monocle_gene_annotation$gene_short_name
  monocle_pd <- new("AnnotatedDataFrame", data = monocle_sample_sheet)
  monocle_fd <- new("AnnotatedDataFrame", data = monocle_gene_annotation)
  # monocle_data <- newCellDataSet(as(expr_matrix, "sparseMatrix"),
  monocle_data <- newCellDataSet(expr_matrix,
                                 phenoData = monocle_pd,
                                 featureData = monocle_fd,
                                 lowerDetectionLimit=1,
                                 expressionFamily=negbinomial.size()) # "Slightly less accurate for differential expression than negbinomial(), but much, much faster." 
  monocle_data <- estimateSizeFactors(monocle_data)
  # NEED TO SUPPRESS STDERR MESSAGES - estimateDispersons popped up a message
  monocle_data <- estimateDispersions(monocle_data)
  monocle_expressed_genes <- row.names(monocle_data)
  monocle_disp_table <- dispersionTable(monocle_data)
  monocle_ordering_genes <- subset(monocle_disp_table,
                                   mean_expression >= 0.1 &
                                   dispersion_empirical >= 1 * dispersion_fit)$gene_id
  monocle_ordering_genes <- monocle_disp_table$gene_id
  monocle_data <- setOrderingFilter(monocle_data, monocle_ordering_genes)
  monocle_data <- reduceDimension(monocle_data, max_components = 2)
  monocle_data <- orderCells(monocle_data, reverse=FALSE)
  # Return
  ##monocle_data
#}

p <- plot_cell_trajectory(monocle_data, color_by=NULL, show_tree=T, show_backbone=T) +
  theme_grey(base_size=20)

# NOTE: Must specify as .png for output to be interpreted
plot_tmp_filename <- "tempggplot.png"
ggsave(filename=plot_tmp_filename, plot=p, width=7, height=5, dpi=100)

p_raw <- readBin(plot_tmp_filename, 'raw', file.info(plot_tmp_filename)$size)
plot_encoded <- base64enc::base64encode(p_raw)
output <- list(
 # name = value,
 # results = list 
 results = list(plot=unbox(plot_encoded))
)

# End STDOUT suppression
sink()
#class(plot_encoded)
#plot_encoded

# Write out JSON with "results"::"plot"
write_json(output, stdout(), simplifyVector=T)
