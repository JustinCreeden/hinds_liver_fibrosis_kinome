library(KEGGREST)
library(KEGGgraph)
pathways = c("mapk" = "hsa04010", "pi3k" = "hsa04151", "ras" = "hsa04014", "egfr" = "hsa01521")

pathway_data = purrr::map(pathways, function(.x){
  
})