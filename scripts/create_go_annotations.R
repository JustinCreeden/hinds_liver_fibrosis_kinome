source(here::here("scripts", "functions.R"))
library(org.Hs.eg.db)
library(dplyr)
go_all = create_go_annotation(org.Hs.eg.db)

saveRDS(go_all, here::here("data", "outputs", "rds_files", "go_all.rds"))

go_bp = create_go_annotation(org.Hs.eg.db, "BP")
saveRDS(go_bp, here::here("data", "outputs", "rds_files", "go_bp.rds"))

go_mf = create_go_annotation(org.Hs.eg.db, "MF")
saveRDS(go_mf, here::here("data", "outputs", "rds_files", "go_mf.rds"))

go_cc = create_go_annotation(org.Hs.eg.db, "CC")
saveRDS(go_cc, here::here("data", "outputs", "rds_files", "go_cc.rds"))

library(reactome.db)
all_reactome = keys(reactome.db)
reactome_df = AnnotationDbi::select(reactome.db, keys = all_reactome,
                     columns = c("PATHNAME", "REACTOMEID")) %>%
  dplyr::filter(grepl("HSA", REACTOMEID))
reactome_split = split(reactome_df$ENTREZID, reactome_df$REACTOMEID)
reactome_split = purrr::map(reactome_split, unique)
reactome_id = unique(reactome_df[, c("PATHNAME", "REACTOMEID")])
reactome_id2 = reactome_id$PATHNAME
reactome_id2 = gsub("Homo sapiens: ", "", reactome_id2)
names(reactome_id2) = reactome_id$REACTOMEID

reactome_annotation = categoryCompare2::annotation(annotation_features = reactome_split,
                                 annotation_type = "reactome",
                                 description = reactome_id2[names(reactome_split)],
                                 feature_type = "ENTREZID")
saveRDS(reactome_annotation, here::here("data", "outputs", "rds_files", "reactome_annotation.rds"))
