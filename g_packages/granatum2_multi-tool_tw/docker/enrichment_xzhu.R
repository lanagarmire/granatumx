# General processing libraries
library(tidyr)
library(dplyr)
library(tibble)
library(stringr)
library(purrr)
library(forcats)
library(ggplot2)

# Specific libraries for this module
library(fgsea)
library(KEGG.db)
library(GO.db)
library(org.Hs.eg.db)
library(org.Mm.eg.db)

source('./granatum_sdk.R')


# Parameters:
# genes = list of DE genes
# scores = z scores from DE analysis
# species = mouse | human
# title = e.g., "cluster_2 vs. all other: Z"

do_kegg_analysis <- function(genes_, scores_, species_, title_, temp_plot_file_name = "temp_kegg_plot.png") {
  species_code <- switch(species_, mouse = 'Mm', human = 'Hs')
  species_code_long <-
    switch(species_, mouse = 'mmu', human = 'hsa')

  kegg_pathways <- mappedkeys(KEGGPATHID2EXTID) %>% str_subset(species_code_long) %>% str_replace(sprintf('^%s(.*)$', species_code_long), '\\1')
  kegg_names <- KEGGPATHID2NAME[kegg_pathways] %>% as.list %>% simplify
  kegg_to_eg <- KEGGPATHID2EXTID[kegg_pathways %>% {sprintf('%s%s', species_code_long, .)}] %>% as.list
  names(kegg_to_eg) <- kegg_names

  SYMBOL2EG <- eval(parse(text = sprintf('org.%s.egSYMBOL2EG', species_code)))
  GO2ALLEGS <- eval(parse(text = sprintf('org.%s.egGO2ALLEGS', species_code)))

  names(scores_) <- genes_
  genes <- intersect(genes_, mappedkeys(SYMBOL2EG))
  scores_ <- scores_[genes]
  gene_entrez <- genes %>% SYMBOL2EG[.] %>% as.list %>% map(~ .[1]) %>% simplify
  names(scores_) <- gene_entrez
  dput(genes)
  fgseaRes <- fgsea(kegg_to_eg, scores_, nperm = 10000)
  full_dat_set <- fgseaRes %>% as.data.frame %>% as_data_frame %>% arrange(-abs(NES)) %>% mutate(pathway = fct_inorder(pathway))
  #full_dat_set <- fgseaRes %>% as.data.frame %>% as_data_frame %>% arrange(-abs(ES)) %>% mutate(pathway = fct_inorder(pathway))
  ggdat <- full_dat_set %>% head(20)

  #print("DEBUG KEGG")
  #print(head(full_dat_set))
  #print(head(full_dat_set$NES))

  # Generate and return raw plot data
  p <- ggplot(ggdat) +
    geom_point(aes(x = pathway,
                   y = abs(NES),
                   #y = abs(ES),
                   size = size)) +
    #labs(title = paste("KEGG:", title_),
    labs(title = paste("KEGG ", title_),
         x = 'Pathway',
         y = 'Absolute Normalized Enrichment Score') +
    scale_size_continuous(name = 'Pathway\nsize') +

    theme_grey(base_size = 12) +
    theme(axis.text.x = element_text(angle = -30, hjust = 0),
          plot.margin = margin(l = 10, r = 10, t = 10, b = 10))

  out_tab <- data.frame(Identifier=full_dat_set$pathway, Score=abs(full_dat_set$NES))
  #out_tab <- data.frame(Identifier=full_dat_set$pathway, Score=abs(full_dat_set$ES))

  # Return
  list(p=p, enrichmentTable=out_tab)
}

do_go_analysis <- function(genes_, scores_, species_, title_, temp_plot_file_name = "temp_go_plot.png") {
  species_code <- switch(species_, mouse = 'Mm', human = 'Hs')

  SYMBOL2EG <- eval(parse(text = sprintf('org.%s.egSYMBOL2EG', species_code)))
  GO2ALLEGS <- eval(parse(text = sprintf('org.%s.egGO2ALLEGS', species_code)))

  names(scores_) <- genes_
  genes <- intersect(genes_, mappedkeys(SYMBOL2EG))

  scores_ <- scores_[genes]

  gene_entrez <-
    genes %>% SYMBOL2EG[.] %>% as.list %>% map(~ .[1]) %>% simplify

  names(scores_) <- gene_entrez

  go_terms <-
    intersect(mappedkeys(GO2ALLEGS), mappedkeys(GOTERM))
  go_to_eg <- GO2ALLEGS[go_terms] %>% as.list
  names(go_to_eg) <-
    names(go_to_eg) %>% GOTERM[.] %>% as.list %>% map(~ .@Term) %>% simplify
  fgseaRes <-
    fgsea(
      go_to_eg,
      scores_,
      nperm = 10000,
      minSize = 50,
      maxSize = 500
    )

  full_dat_set <-
    # NES was causing "Inf" values to be calculated for a case when using GO;
    # but ES results do not look good so keeping with NES;
    # it also looks like if too few genes then could throw an error - the marker
    # selection step should output a full list of genes plus scores [20180717 tkwolf]
    fgseaRes %>% as.data.frame %>% as_data_frame %>% arrange(-abs(NES)) %>% mutate(pathway =
    #fgseaRes %>% as.data.frame %>% as_data_frame %>% arrange(-abs(ES)) %>% mutate(pathway =
                                                                                     fct_inorder(pathway))
  #print("DEBUG GO")
  #print(head(full_dat_set))
  #print(head(full_dat_set$NES))
  #print(head(full_dat_set$ES))

  ggdat <- full_dat_set %>% head(20)

  # Generate and return raw plot data
  p <- ggplot(ggdat) +
    geom_point(aes(x = pathway,
                   y = abs(NES),
                   #y = abs(ES),
                   size = size)) +
    #labs(title = paste("GO:", title_),
    labs(title = paste("GO ", title_),
         x = 'Gene set',
         y = 'Absolute Normalized Enrichment Score') +
    scale_size_continuous(name = 'Gene set\nsize') +

    theme_grey(base_size = 12) +
    theme(axis.text.x = element_text(angle = -30, hjust = 0),
          plot.margin = margin(l = 10, r = 10, t = 10, b = 10))

  out_tab <- data.frame(Identifier=ggdat$pathway, Score=abs(ggdat$NES))
  #out_tab <- data.frame(Identifier=ggdat$pathway, Score=abs(ggdat$ES))

  # Return
  list(p=p, enrichmentTable=out_tab)
}


#set.seed(jsonInput$seed)

# Method: KEGG | GO
enrichmentMethod <- gn_get_arg('geneSetDatabase')

# WARNING! Expects the first element in the list to have the gene/value
# records that will be used for enrichment analysis
geneScores <- gn_get_import('genesAndScores')
# We would need another way to get the label for this data, e.g., Cluster 2 vs. rest. [20180717 tkwolf]
#title <- names(geneScores)
title <- '' # KEGG or GO will be prefixed to this
species <- gn_get_arg('species')

# No exports for now until we decide on a way to store it since choices at
# the moment are assay, geneMeta, and sampleMeta; we might consider using
# geneMeta by assigning the pathway(s)/term(s) associated with each gene if
# we want to use the pathways/terms in another processing step or to create
# a something like a setMeta type.
# finalOutput[["exports"]] <- list()

p <- NULL
enrichmentTable <- NULL
plotDescription <- NULL
# TEST
# enrichmentMethod = "GO"
if (enrichmentMethod == "KEGG") {
  analysis <-
    do_kegg_analysis(names(geneScores),
                     unlist(geneScores),
                     species,
                     title,
                     "temp_kegg_plot.png")
  p <- analysis$p
  enrichmentTable <- analysis$enrichmentTable
  plotDescription <- "Plot of enriched KEGG pathways given differential expression results."
} else if (enrichmentMethod == "GO") {
  analysis <-
    do_go_analysis(names(geneScores),
                   unlist(geneScores),
                   species,
                   title,
                   "temp_go_plot.png")
  p <- analysis$p
  enrichmentTable <- analysis$enrichmentTable
  plotDescription <- "Plot of enriched GO terms given differential expression results."
} else {
  # Error condition
}

gn_add_ggplot_to_results(p, 'Enrichment plot', dpi=150, height=400)
gn_add_dataframe_to_results(enrichmentTable, 'Enrichment table', rowname_header='Rank')

gn_export_dynamically(enrichmentTable, 'Enrichment table')

gn_commit()
