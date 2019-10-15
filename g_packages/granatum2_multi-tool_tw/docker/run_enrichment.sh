#!/usr/bin/bash
# Note that the shebang above is necessary to run in the Docker container
# Input: JSON including DE results and enrichment method, e.g., KEGG
# Output: JSON with table and plot of enriched pathways/terms
Rscript enrichment.R #2>/dev/null
