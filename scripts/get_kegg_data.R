library(KEGGREST)
source(here::here("scripts", "functions.R"))
hsa_pathways <- keggLink("hsa", "pathway")
names(hsa_pathways) <- substring(names(hsa_pathways), 6)
hsa_pathways <- substring(hsa_pathways, 5)
kegg_annotation <- list(annotation = split(hsa_pathways, names(hsa_pathways)))
kegg_annotation$annotation <- lapply(kegg_annotation$annotation, function(x){names(x) <- NULL; x})
 
kegg_paths <- names(kegg_annotation$annotation)
kegg_desc <- keggList("pathway", "hsa")
names(kegg_desc) <- substring(names(kegg_desc), 6)
kegg_annotation$description <- kegg_desc[names(kegg_annotation$annotation)]

kegg_annotation2 = create_kegg_annotation(kegg_annotation)
saveRDS(kegg_annotation2, here::here("data", "outputs", "rds_files", "kegg_annotation.rds"))
# last run on 2021-10-27