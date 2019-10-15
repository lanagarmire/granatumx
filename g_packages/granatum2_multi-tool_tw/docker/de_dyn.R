# 20180302 twolf
# Load dependencies
suppressPackageStartupMessages({
  library(jsonlite)
})

in_json <- read_json(
  'stdin',
  simplifyVector = T
)

# Input number of groups to calculate the number of outputs
# Each group is compared to all other samples (all other groups)
groups <- in_json$groups
# Need actual samples in matrix, since will only use groups from those
expr_samples <- in_json$assay$sampleIds
# Need method to describe output
#   NODES, limma, edgeR, SCDE
selected_de <- in_json$deMethod

# Getting groups for input samples
expr_samples <- expr_samples[expr_samples %in% names(groups)]
expr_groups <- NULL
for (s in expr_samples) {
  expr_groups <- c(expr_groups, groups[s])
}
expr_groups <- as.vector(unlist(expr_groups))

# If two groups then format:
#   sprintf('%s vs. %s', pair[1], pair[2])
# If >2 groups then format:
#   sprintf('%s vs. all other', l)
# Modifying the function that is used in the de.R script
#   to keep things consistent.
get_pairs <- function (vec) {
  vec <- factor(vec)
  lvls <- levels(vec)

  # Consider revising below - how should we pass this message back?
  if (length(lvls) < 2)
    stop('Need at least two groups for DE.')
  
  pair_tags <- NULL
  
  # if (length(lvls) == 2 || pairwise) {
  if (length(lvls) == 2) {
    pairs <- as.list(as.data.frame(combn(lvls, 2)))
    for (pair in pairs) {
      # Already checked:
      # vec2 <- vec[vec %in% pair]
      # mat2 <- mat[, vec %in% pair]
      pair_tags <- sprintf('%s vs. %s', pair[1], pair[2])
    }
  } else {
    for (l in lvls) {
      lvls2        <- rep(2, length(lvls))
      names(lvls2) <- lvls
      lvls2[l]     <- 1
      vec2         <- lvls2[vec]
      
      pair_tags <- c(pair_tags, sprintf('%s vs. all other', l))
    }
  }
  # Return
  pair_tags
}

# For types: https://developer.mozilla.org/en-US/docs/Web/HTML/Element
output <- list()
get_desc <- function (method) {
  # Simple for now, later add column meanings per method
  desc <- switch(
    method,
    'NODES' = "Generated from NODES.",
    'SCDE' = "Generated from SCDE.",
    'edgeR' = "Generated from edgeR.",
    'limma' = "Generated from limma."
  )
  desc
}
get_fields <- function (method) {
  # Careful to keep the order the same as JSON output that main DE will output
  # Using the JSON output order from the main DE script
  # NODES e.g.:
  #   {"results":{"cluster_0 vs. all other":{"data":[{"Fisher":3.5106e-09,"qvalues":1.6161e-06,"Z":-3915.6035,"abs_Z":3915.6035,"Gene":"ACTB"}, ...
  #   becomes
  #   c("Fisher","qvalues","Z","abs_Z")
  desc <- switch(
    method,
    'NODES' = c("Fisher","qvalues","Z","abs_Z"),
    'SCDE' = c("lb","mle","ub","ce","Z","cZ","p.values","p.values.adj"),
    'edgeR' = c("logFC","logCPM","PValue","FDR","Z"),
    'limma' = c("logFC","AveExpr","t","P.Value","adj.P.Val","B","Z")
  )
  desc
}
all_fields <- get_fields(selected_de)
all_pairs <- get_pairs(expr_groups)
for (i in 1:length(all_pairs)) {
  t <- all_pairs[i]
  output$results[[i]] <- list(
    type = unbox('table'),
    label = unbox(t),
    extractFrom = unbox(t),
    description = unbox(get_desc(selected_de))
  )
}
exports_index <- 1
for (i in 1:length(all_pairs)) {
  p <- all_pairs[i]
  for (j in 1:length(all_fields)) {
    f <- all_fields[j]
    pair_field = sprintf("%s: %s",p,f)
    output$exports[[exports_index]] <- list(
      kind = unbox('geneMeta'),
      extractFrom = unbox(pair_field),
      name = unbox(pair_field)
    )
    exports_index <- exports_index + 1
  }
}

# jsonlite::toJSON(output)
write_json(output, stdout(), simplifyVector=T)
