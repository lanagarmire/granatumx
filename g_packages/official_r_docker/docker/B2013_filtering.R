plot_tmp_filename <- '/tmp/plot.png'

suppressMessages(library(jsonlite))
suppressMessages(library(reshape2))
suppressMessages(library(tidyverse))

input <- read_json('stdin', simplifyVector=T)

assay <- input$assay

mat <- assay$matrix

xx <- input$expressedAbove
yy <- input$inPercentage/100

message(paste(capture.output(nrow(mat)), collapse='\n'))
message(paste(capture.output(ncol(mat)), collapse='\n'))
message(paste(capture.output(xx), collapse='\n'))
message(paste(capture.output(yy), collapse='\n'))

keptGenes <- apply(mat, 1, function(x)mean(x > xx) > yy)



# --------- CV2 plot ---------

tb <- melt(mat, varnames=c('gene', 'sample'), value.name='value')

avg <- apply(mat, 1, mean)
cv2 <- apply(mat, 1, var) / avg^2

p <- ggplot() +
  geom_point(aes(avg, cv2, color=keptGenes)) +
  scale_color_discrete(name='', labels=c('Removed genes', 'Kept genes')) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x='Average expression across all cells', y='Variance / (Avg. ^ 2)')

ggsave(plot_tmp_filename, p, width=7.5, height=5, dpi=100)

p_raw <- readBin(plot_tmp_filename, 'raw', file.info(plot_tmp_filename)$size)

plot_encoded <- base64enc::base64encode(p_raw)

# ----------------------------


assay$matrix <- mat[keptGenes,]
assay$geneIds <- assay$geneIds[keptGenes]


output <- list(
  exports = list(
    filteredAssay = assay
  ),
  results = list(
    numGenesBefore=unbox(nrow(mat)),
    numGenesAfter=unbox(nrow(assay$matrix)),
    plot=plot_encoded
  )
)

write_json(output, stdout())
