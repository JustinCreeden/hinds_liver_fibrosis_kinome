# Run script to completely recreate R analyses from the manuscript.
# {future} package is used to run each script in it's own environment.
# If you want to see the records of what was run or errors from a 
#   particular script running, then you should use `value`.
#   For example, for the first script to run:
#   kea3_enrichment = future(source(here::here("scripts", "kea3_runs.R"))) # runs it
#   value(kea3_enrichment) # shows any output from running it
#   
# Everything that was previously generated is included in the output
# and reports directories,
# so these are only necessary if you want to update something.
# See the README.md file for documentation of what is in each directory.
library(future)

# download and unzip the HPA consensus data
download.file("https://www.proteinatlas.org/download/rna_tissue_consensus.tsv.zip",
              destfile = here::here("data", "inputs", "hpa", "rna_tissue_consensus.tsv.zip"))
unzip(here::here("data", "inputs", "hpa", "rna_tissue_consensus.tsv.zip"),
      exdir = here::here("data", "inputs", "hpa"))

kea3_enrichment = future(source(here::here("scripts", "kea3_runs.R")))
kea3_networks = future(source(here::here("scripts", "create_kea3_network.R")))
go_annotation = future(source(here::here("scripts", "create_go_annotations.R")))
kegg_data = future(source(here::here("scripts", "get_kegg_data.R")))
kegg_graphs = future(source(here::here("scripts", "pull_kegg_pathways.R")))

kegg_network_graph = future(rmarkdown::render(here::here("reports", "kegg_pathway_graphing.Rmd"), output_format = "word_document"))
pvalue_adjust = future(rmarkdown::render(here::here("reports", "pvalue_adjustments.Rmd")))
kea3_enetwork = future(rmarkdown::render(here::here("reports", "kea3-enrichment_network.Rmd")))
pathway_enrichment = future(rmarkdown::render(here::here("reports", "pathway_enrichment.Rmd")))
fib_normal_raw_correlation = future(rmarkdown::render(here::here("reports", "fib_normal_raw_correlation.Rmd")))
fib_normal_signal_correlation = future(rmarkdown::render(here::here("reports", "fib_normal_signal_correlation.Rmd")))
hpa_consensus = future(rmarkdown::render(here::here("reports", "exploring_hpa_consensus.Rmd")))
