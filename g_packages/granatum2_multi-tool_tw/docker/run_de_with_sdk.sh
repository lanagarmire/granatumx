#!/usr/bin/bash
# Note that the shebang above is necessary to run in the Docker container
# Run differential expression
# Pass JSON to R file that loads functions in differential_expression/*.R files,
# parses the JSON input, runs the processing, and outputs results in JSON.
#
#Rscript de.R 2>/dev/null
Rscript de_with_sdk.R
