# Run script to completely recreate R analyses from the manuscript.
# {future} package is used to run each script in it's own environment.
# If you want to see the records of what was run or errors from a 
#   particular script running, then you should examine it's output.
#   For example, for the first script to run:
#   kea3_enrichment = r(function() source(here::here("scripts", "kea3_runs.R"))) # runs it
#   kea3_enrichment # shows any output from running it
#   ..., show = TRUE) # will show what is going on in the r sub-process
# Everything that was previously generated is included in the output
# and reports directories,
# so these are only necessary if you want to update something.
# See the README.md file for documentation of what is in each directory.
library(callr)

# download and unzip the HPA consensus data
download.file("https://www.proteinatlas.org/download/rna_tissue_consensus.tsv.zip",
              destfile = here::here("data", "inputs", "hpa", "rna_tissue_consensus.tsv.zip"))
unzip(here::here("data", "inputs", "hpa", "rna_tissue_consensus.tsv.zip"),
      exdir = here::here("data", "inputs", "hpa"))

kea3_enrichment = r(function() source(here::here("scripts", "kea3_runs.R")))
kea3_networks = r(function() source(here::here("scripts", "create_kea3_network.R")))
go_annotation = r(function() source(here::here("scripts", "create_go_annotations.R")))
kegg_data = r(function() source(here::here("scripts", "get_kegg_data.R")))
kegg_graphs = r(function() source(here::here("scripts", "pull_kegg_pathways.R")))

kegg_network_graph = r(function() rmarkdown::render(here::here("reports", "kegg_pathway_graphing.Rmd"), output_format = "word_document"), show = TRUE)
pvalue_adjust = r(function() rmarkdown::render(here::here("reports", "pvalue_adjustment.Rmd")), show = TRUE)
kea3_enetwork = r(function() rmarkdown::render(here::here("reports", "kea3-enrichment_network.Rmd")), show = TRUE)
pathway_enrichment = r(function() rmarkdown::render(here::here("reports", "pathway_enrichment.Rmd")), show = TRUE)
fib_normal_raw_correlation = r(function() rmarkdown::render(here::here("reports", "fib_normal_raw_correlation.Rmd")), show = TRUE)
fib_normal_signal_correlation = r(function() rmarkdown::render(here::here("reports", "fib_normal_signal_correlation.Rmd")), show = TRUE)
hpa_consensus = r(function() rmarkdown::render(here::here("reports", "exploring_hpa_consensus.Rmd")), show = TRUE)
