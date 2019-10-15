# Load dependencies
library(tidyr)
library(dplyr)
library(tibble)
library(MetaDE)
library(NODES)
library(limma)
library(edgeR)
library(scde)

source('./granatum_sdk.R') # Uses GranatumX SDK

# Start of DE functions definitions
## NODES
#' Nonparametric differential expression analysis for scRNA-seq data.
#' @param data \code{pQ} normalized single cell data.
#' @param group An array of strings containing unique group identifiers for cells. For example \code{c('A','A','B','B')}.
#' @param r A numeric value indicating number of permutations of outcomes to be done for
#' generting the empirical distribution for \code{D} statistic. Default value is \code{20}.
#' @param smooth_points A numeric value indicating number of bins for smoothing, default is \code{10000}.
#' @param zper Indicating quantile of pooled standard errors; default is \code{0.5}.
#' @return A table reporting \code{p}-values and \code{q}-values with original ordering of genes retained.

#' data(data_Trapnell) # load the Trapnell data
#' norm_Data <- pQ(data_Trapnell) # pQ normalize the data
#' grp <- c(rep('T0',75),rep('T24',71)) # Group assignment
#' Res <- NODES(norm_Data[1:10,],grp) # Compute p-value just for 10 genes

NODES <- function(data, group, r = 20, smooth_points = 10000, zper = 0.5) {
  # Put into dependencies above
  # # install require libraries
  # if (!require("MetaDE"))
  #   install.packages("MetaDE", repos = "http://cran.us.r-project.org")
  # 
  # require(MetaDE)

  # This part is for identifying the groups
  indices <- list()

  U <- unique(group)

  for (i in 1:length(U)) {
    indices[[as.character(U[i])]] <- grep(U[i], group)
  }

  # length of group 1
  n1 <- length(indices[[U[1]]])
  n2 <- length(indices[[U[2]]])

  # Getting Noise distribution from NOISeq and estimate p values

  Zr <- NULL
  for (i in 1:r) {

    # print(paste("Randomization run =", i))

    mipermu <- sample(1:(n1 + n2))  ## randomize labels

    mipermu <- data[, mipermu]  ## randomize matrix columns accordingly

    mean1 <- rowMeans(mipermu[, 1:n1])  ## get the means for random group 1
    mean2 <- rowMeans(mipermu[, (n1 + 1):(n1 + n2)])  ## get the means for random group 2

    sd1 <- apply(mipermu[, 1:n1], 1, sd)  ## sd for group 1
    sd2 <- apply(mipermu[, (n1 + 1):(n1 + n2)], 1, sd)  ## sd for group 2

    myparam <- list(n = c(n1, n2), sd = cbind(sd1, sd2))

    MDperm <- MDbio(dat = cbind(mean1, mean2), param = myparam, a0per = zper)

    Zr <- cbind(Zr, MDperm$D)

  }

  # Estimating noise density using Gaussian kernel
  cat("\nSmoothing the noise density...\n")
  dF <- approxfun(density(as.vector(Zr), n = smooth_points))

  # Getting stat for all genes

  mean1 <- rowMeans(as.matrix(data[, indices[[U[1]]]]))
  mean2 <- rowMeans(as.matrix(data[, indices[[U[2]]]]))

  sd1 <- apply(as.matrix(data[, indices[[U[1]]]]), 1, sd)
  sd2 <- apply(as.matrix(data[, indices[[U[2]]]]), 1, sd)

  myparam <- list(n = c(n1, n2), sd = cbind(sd1, sd2))

  Ds <- MDbio(dat = cbind(mean1, mean2), param = myparam, a0per = zper)

  Zs <- Ds$D
  Zrs <- Ds$Dr


  # Estimating p values only from noise distribution

  prob_1 <- apply(as.matrix(Zs), 1, function(x) den(dF, as.vector(Zr),
    x, smooth_points))



  # getting wilcoxon p values
  cat("\nComputing Wilcoxon p values...\n")
  prob_2 <- apply(data, 1, function(x) wilcox.test(x[indices[[U[1]]]],
    x[indices[[U[2]]]])$p.value)

  # Fisher's method to combine p values
  cat("\nCombining p values using Fisher's method...\n")

  together <- list()


  together[["p"]] <- as.matrix(cbind(prob_1, prob_2))
  META <- MetaDE::MetaDE.pvalue(together, meta.method = "Fisher")
  PVAL <- META$meta.analysis$pval
  # names(PVAL)<-rownames(data)

  # fdr
  fdr <- p.adjust(PVAL, method = "fdr")

  # prepare final result
  res <- data.frame(cbind(pvalues = PVAL, qvalues = fdr, Z = Zrs))

  rownames(res) <- rownames(data)

  res$abs_Z <- abs(res$Z)

  res <- res[order(-abs(res$Z)), ]

  cat("\nCompleted successfully.\n")

  # return
  return(res)
}

# Tail probability estimation
den <- function(approx, obs, val, points) {
  # if(val>=max(obs) || val <= min(obs))
  if (val >= max(obs)) {
    pF <- 1/(1 + points)
  } else {
    # pF <-
    # min(integrate(approx,val,max(obs))$value,integrate(approx,min(obs),val)$value)
    pF <- integrate(approx, val, max(obs), subdivisions=2000)$value
  }

  return(pF)  # 2*pF otherwise
}

## This code is a shadow of NOISeqbio
MDbio <- function(dat = dat, param = NULL, a0per = 0.5) {
  ddr <- (dat[,1]-dat[,2]) # this is for two tailed.
  dd <- abs(dat[, 1] - dat[, 2])


  # sd.D = sqrt(param$sd[,1]^2/sqrt(param$n[1]) +
  # param$sd[,2]^2/sqrt(param$n[2]))
  sd.D <- sqrt((param$sd[, 1]^2/param$n[1]) + (param$sd[, 2]^2/param$n[2]))

  a0per <- as.numeric(a0per)
  a0.D <- quantile(sd.D, probs = a0per, na.rm = TRUE)

  dd <- dd/(a0.D + sd.D)

  # Results
  list(D = dd, Dr=ddr)
}


## limma
do_limma_for_two_groups <- function (mat, vec, n_cores)  {
  vec <- paste('group_', as.character(vec), sep='')
  design <- model.matrix(~0+vec)
  lvls <- vec %>% factor %>% levels
  colnames(design) <- lvls
  contrast <- makeContrasts(contrasts=sprintf('%s - %s', lvls[1], lvls[2]), levels=design)
  fit <- lmFit(mat, design)
  fit <- contrasts.fit(fit, contrast)
  fit <- eBayes(fit)
  res <- topTable(fit, number=Inf)
  res$Z <- res$logFC

  res
}


## EdgeR
do_edgeR_for_two_groups <- function (mat, vec, n_cores)  {
  cds <- DGEList(mat, group=vec)
  cds <- cds[rowSums(1e+06 * cds$counts/expandAsMatrix(cds$samples$lib.size, dim(cds)) > 1) >= 3, ]
  cds <- calcNormFactors( cds )
                                        #plotMDS( cds , main = "MDS Plot for Count Data", labels = sub("^.*(.)..","\\1",colnames( cds$counts )) )
  de <- exactTest(cds, dispersion=0.2)

  res <- topTags(de, Inf)$table

  res$Z <- res$logFC

  res
}

do_nodes_for_two_groups <- function (mat, vec, n_cores)  {
  NODES(mat, as.character(vec))
}


## SCDE
do_scde_for_two_groups <- function (mat, vec, n_cores)  {
  scde_counts             <- data.frame(mat)
  scde_sg                 <- factor(vec)

  names(scde_sg)          <- colnames(scde_counts)

  scde_cd                 <- scde_counts
  scde_o_ifm              <- scde.error.models(counts = scde_cd,
                                               groups = scde_sg,
                                               n.cores = n_cores,
                                               min.size.entries = min(2000, nrow(mat)),
                                               threshold.segmentation = TRUE,
                                               save.model.plots = FALSE,
                                               verbose = 1)
  scde_valid_cells        <- scde_o_ifm$corr.a > 0
  scde_o_ifm              <- scde_o_ifm[scde_valid_cells,]
  scde_o_prior            <- scde.expression.prior(models = scde_o_ifm, counts = scde_cd, length.out = 400, show.plot = FALSE)
  scde_sg_valid           <- scde_sg[rownames(scde_o_ifm)]

  scde_ediff              <- scde.expression.difference(scde_o_ifm,
                                                        scde_cd, scde_o_prior,
                                                        groups = scde_sg_valid,
                                                        n.randomizations = 100,
                                                        n.cores = n_cores,
                                                        verbose = 0)

  p.values <- 2*pnorm(abs(scde_ediff$Z),lower.tail=F) # 2-tailed p-value
  p.values.adj <- 2*pnorm(abs(scde_ediff$cZ),lower.tail=F) # Adjusted to control for FDR

  scde_ediff %>% as.data.frame %>% rownames_to_column('Gene') %>%
    mutate(p.values = p.values) %>%
    mutate(p.values.adj = p.values.adj) %>%
    arrange(p.values.adj) %>%
    as.data.frame %>%
    column_to_rownames('Gene')
}


## Function to run the selected method function
do_diff_exp <- function (mat, vec, pairwise=F, n_cores, method)  {
  vec <- factor(vec)
  lvls <- levels(vec)

  # Consider revising below - how should we pass this message back?
  if (length(lvls) < 2) stop('Need at least two groups for DE.')

  output <- list()

  if (length(lvls) == 2 || pairwise) {
    pairs <- as.list(as.data.frame(combn(lvls, 2)))
    for (pair in pairs) {
      vec2 <- vec[vec %in% pair]
      mat2 <- mat[,vec %in% pair]

      output[[sprintf('%s vs. %s', pair[1], pair[2])]] <-
        switch(method,
               'NODES' = do_nodes_for_two_groups(mat2, vec2, n_cores),
               'SCDE' = do_scde_for_two_groups(mat2, vec2, n_cores),
               'edgeR' = do_edgeR_for_two_groups(mat2, vec2, n_cores),
               'limma' = do_limma_for_two_groups(mat2, vec2, n_cores))
    }
  } else {
    for (l in lvls) {
      lvls2        <- rep(2, length(lvls))
      names(lvls2) <- lvls
      lvls2[l]     <- 1
      vec2         <- lvls2[vec]

      output[[sprintf('%s vs. all other', l)]] <-
        switch(method,
               'NODES' = do_nodes_for_two_groups(mat, vec2, n_cores),
               'SCDE' = do_scde_for_two_groups(mat, vec2, n_cores),
               'edgeR' = do_edgeR_for_two_groups(mat, vec2, n_cores),
               'limma' = do_limma_for_two_groups(mat, vec2, n_cores))
    }
  }

  output
}


# Input parameters/data and calculate DE table
# Import assay and selected option(s)
# TEST
# in_json <- read_json('input_with_genes_de.json', simplifyVector=T)
# in_json <- read_json('stdin', simplifyVector=T)

# Set the seed for reproducibility
rand_seed <- gn_get_arg('seed')
set.seed(rand_seed)
# Tested seed (96813) using diff (temp1 vs. temp2, etc.) of output using 6k subset of features by the method of...
#   NODES: OK
#   limma: OK
#   edgeR: OK
#   SCDE: OK

# Method: NODES, limma, edgeR, SCDE
selected_de <- gn_get_arg('deMethod')

# Number of compute cores to use
selected_cores <- gn_get_arg('deCores')

# Normalized assay
assay <- gn_get_import('assay')
expr_matrix <- assay$matrix
# Sample IDs in matrix; assumes the order is same as order of columns in matrix
expr_samples <- assay$sampleIds
# Gene IDs; needed for output
g_ids <- assay$geneIds

# Groups to compare, e.g., from clustering cells; only use these samples
sampleMeta <- gn_get_import('sampleMeta')
#groups <- sampleMeta$groups
groups <- sampleMeta

expr_matrix <- expr_matrix[,expr_samples %in% names(groups)]
expr_samples <- expr_samples[expr_samples %in% names(groups)]
expr_groups <- NULL
for (s in expr_samples) {
  expr_groups <- c(expr_groups, groups[s])
}
expr_groups <- as.vector(unlist(expr_groups))

# TEST
# expr_matrix <- head(expr_matrix, 2000)
# g_ids <- g_ids[1:2000]
# selected_de <- 'NODES'
# # selected_de <- 'limma'
# # selected_de <- 'edgeR'
# # selected_de <- 'SCDE'

rownames(expr_matrix) <- g_ids

# SCDE requires integer values - expects raw counts (not, for example, FPKM).
# As a temporary fudge for SCDE we are recasting matrix values as integers.
# This should be handled better when we split off each method as their
# own individual Gbox (scripts) with appropriate pre-requisites.
if (selected_de == 'SCDE') {
  expr_matrix <- apply(expr_matrix,2,function(x) {storage.mode(x) <- 'integer'; x})
}

out_tab <- do_diff_exp(expr_matrix, expr_groups, pairwise=F, n_cores=selected_cores, selected_de)

# Output allowing for dynamic UI in these results output (revision 3)
# Results: (displayed results for the "second page")
#   * type : png | text | iframe | scatterplot | table
#   * data : (huge data)
#   * label : '....'
#   --- below are depending on type ---
#   * width
#   * height
#   * description
output = list()
# Results section
# output[["results"]] = list()
for (i in 1:length(names(out_tab))) {
  output[["results"]][[i]] <- out_tab[i]
  names(output[["results"]][[i]]) <- "data"
  output[["results"]][[i]][["data"]][["Gene"]] <- rownames(out_tab[[i]]) # Careful that this is different than SDK-defined rowname_header [20180720 tkwolf]
  rownames(output[["results"]][[i]][["data"]]) <- NULL
  #output[["results"]][[i]][["type"]] <- unbox("table")
  #output[["results"]][[i]][["label"]] <- unbox(names(out_tab)[i])
  #output[["results"]][[i]][["description"]] <- unbox(sprintf("Differential expression results from %s using %s", names(out_tab)[i], selected_de))
  output[["results"]][[i]][["type"]] <- "table"
  output[["results"]][[i]][["label"]] <- names(out_tab)[i]
  output[["results"]][[i]][["description"]] <- sprintf("Differential expression results from %s using %s", names(out_tab)[i], selected_de)
  output[["results"]][[i]][["cols"]] <- list()
  for (j in 1:length(names(out_tab[[i]]))) {
    field_name <- names(out_tab[[i]])[j]
    #output[["results"]][[i]][["cols"]][[j]] <- list("name" = unbox(field_name))
    output[["results"]][[i]][["cols"]][[j]] <- list("name" = field_name)
  }
}
names(output[["results"]]) <- names(out_tab)

# Exports section
output[["exports"]] = list()

# When going to "dynamic" output version, uncomment line below:
# exports_index <- 1

for (p in names(out_tab)) {
  for (f in names(out_tab[[p]])) {

    # When going to "dynamic" output version, uncomment line below:
    # output[["exports"]][[exports_index]] <- list()
    
    pair_field = sprintf("%s: %s",p,f)
    
    # When going to "dynamic" output version, comment or remove block below:
    output[["exports"]][[pair_field]] = as.list(out_tab[[p]][[f]])
    names(output[["exports"]][[pair_field]]) = rownames(out_tab[[p]])
    for (i in (1:length(output[["exports"]][[pair_field]]))) {
      #output[["exports"]][[pair_field]][[i]] <- unbox(output[["exports"]][[pair_field]][[i]])
      output[["exports"]][[pair_field]][[i]] <- output[["exports"]][[pair_field]][[i]]
    }
    
    # When going to "dynamic" output version, uncomment block below below:
    # output[["exports"]][[exports_index]][["data"]] = as.list(out_tab[[p]][[f]])
    # names(output[["exports"]][[exports_index]][["data"]]) = rownames(out_tab[[p]])
    # for (i in (1:length(output[["exports"]][[exports_index]][["data"]]))) {
    #   output[["exports"]][[exports_index]][["data"]][[i]] <- unbox(output[["exports"]][[exports_index]][["data"]][[i]])
    # }
    # output[["exports"]][[exports_index]][["kind"]] <- unbox("geneMeta")
    # output[["exports"]][[exports_index]][["name"]] <- unbox(pair_field)
    # exports_index <- exports_index + 1
    
  }
}

# Output
# TEST
# write_json(output, "temp", simplifyVector=T)
# write_json(output, stdout(), simplifyVector=T)

# Results
for (i in 1:length(names(out_tab))) {
  gn_add_dataframe_to_results(output[["results"]][[i]][["data"]],
                              description=output[["results"]][[i]][["description"]],
                              rowname_header="Record") # Note: this cannot be the same as an existing column name, e.g., "Gene" [20180720 tkwolf]
}

# Exports
for (p in names(out_tab)) {
  for (f in names(out_tab[[p]])) {
    pair_field = sprintf("%s: %s",p,f)
#    for (i in (1:length(output[["exports"]][[pair_field]]))) {
#      gn_export_dynamically(output[["exports"]][[pair_field]][[i]],
#                            sprintf("DE table: %s",pair_field),
#                            kind="geneMeta",
#                            meta=NULL)
#    }
    gn_export_dynamically(output[["exports"]][[pair_field]],
                         sprintf("DE table: %s",pair_field),
                         kind="geneMeta",
                         meta=NULL)
  }
}

gn_commit()
