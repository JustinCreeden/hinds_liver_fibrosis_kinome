library(KEGGREST)
library(KEGGgraph)
library(dplyr)
pathways = c("mapk" = "hsa04010", "pi3k" = "hsa04151", "ras" = "hsa04014", "egfr" = "hsa01521")

pathway_data = purrr::map(pathways, function(.x){
  kgml = keggGet(.x, "kgml")
  graph_df = parseKGML2DataFrame(kgml) %>%
    dplyr::mutate(from = gsub("hsa:", "", from),
                  to = gsub("hsa:", "", to))
  graph_df
})

saveRDS(pathway_data, file = here::here("data", "outputs", "rds_files", "kegg_pathway_graphs.rds"))
