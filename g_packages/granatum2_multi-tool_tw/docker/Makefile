PWD=.
D=docker

# Run commands
R_PT=run_pseudo-time.sh
R_N=run_normalization.sh
# The _DYN should be run first to tell the core what tables to expect
R_DE_DYN=run_de_dyn.sh
R_DE=run_de.sh
R_FG=run_filter_genes.sh
R_EN=run_enrichment.sh

# Usage script should be run by default, but given here for clarification
R_U=run_usage.sh

# Input test JSON files
I_PT=input_test_pseudo-time.json
#
I_N=input_test_normalization.json
#
# I_DE=input_test_de.json
# I_DE_3=input_test_de_3_groups.json
# I_DE_SCDE=input_test_de_scde_small.json
I_DE=input_with_genes_de.json
I_FG=input_with_genes_de_plus.json
I_EN=input_test_enrichment.json

# Using "@" to not echo the command

docker_image: Dockerfile
	# Granatum2 multi-tool by TW
	$(D) --version
	$(D) build -t granatum2_multi-tool_tw $(PWD)

usage:
	$(D) run -i --rm granatum2_multi-tool_tw

test_pseudo-time:
	cat $(I_PT) | $(D) run -i --rm granatum2_multi-tool_tw ./$(R_PT)

test_normalization:
	cat $(I_N) | $(D) run -i --rm granatum2_multi-tool_tw ./$(R_N)

test_de_dyn:
	@cat $(I_DE) | $(D) run -i --rm granatum2_multi-tool_tw ./$(R_DE_DYN)

test_de:
	@cat $(I_DE) | $(D) run -i --rm granatum2_multi-tool_tw ./$(R_DE)

# test_de_3:
# 	cat $(I_DE_3) | $(D) run -i --rm granatum2_multi-tool_tw ./$(R_DE)
# 
# test_de_scde:
# 	cat $(I_DE_SCDE) | $(D) run -i --rm granatum2_multi-tool_tw ./$(R_DE)

test_gene_filtering:
	@cat $(I_FG) | $(D) run -i --rm granatum2_multi-tool_tw ./$(R_FG)

test_enrichment:
	@cat $(I_EN) | $(D) run -i --rm granatum2_multi-tool_tw ./$(R_EN)
