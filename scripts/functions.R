# functions used for reading in and processing data
# for the KEA3 enrichment analysis.

filter_hpa_expression = function(hpa_expression, tissue = NULL, cutoff = NULL, n_sd = 2){
  
  if (!is.null(tissue)) {
    tissue_expression = hpa_expression %>%
      dplyr::filter(Tissue %in% tissue)
  } else {
    tissue_expression = hpa_expression
  }
  tissue_expression = tissue_expression %>%
    dplyr::mutate(logX = log1p(NX))
  
  if (!is.null(cutoff)) {
    tissue_isexpressed = tissue_expression %>%
      dplyr::filter(NX > cutoff)
    return(tissue_isexpressed)
  }
  
  pos_expression = tissue_expression %>%
    dplyr::filter(logX >= 1) %>%
    dplyr::pull(logX)
  
  mode_sd = calculate_mode_sd(pos_expression)
  
  min_val = mode_sd[1] - (n_sd * mode_sd[2])
  
  tissue_isexpressed = tissue_expression %>%
    dplyr::filter(logX >= min_val)
  
  return(tissue_isexpressed)
}

parse_gmt = function(gmt_file){
  gmt = scan(gmt_file, what = "", sep = "\n", quiet = TRUE)
  gmt_split = strsplit(gmt, "\t")
  names(gmt_split) = purrr::map_chr(gmt_split, ~ gsub("_.*", "", .x[[1]]))
  gmt_split = purrr::map(gmt_split, ~ .x[-1])
  gmt_split
}


calculate_mode_sd = function(pos_expression){
  density_expression = stats::density(pos_expression)
  high_y = which.max(density_expression$y)
  mode_val = density_expression$x[high_y]
  
  high_vals = pos_expression[pos_expression >= mode_val]
  sd_high = sqrt(sum((high_vals - mode_val)^2) / (length(high_vals) - 1))
  return(c("mode" = mode_val, "sd" = sd_high))
}

create_kea3_networks = function(){
  kea3_data_dir = here::here("data", "inputs", "kea3")
  all_files = dir(kea3_data_dir, full.names = TRUE, pattern = "gmt$")
  ks_libs = c("Cheng.KSIN", "PTMsigDB", "PhosD.All")
  ppi_libs = c("BioGRID", "MINT", "Mentha", "HIPPIE", "PrePPI", "Cheng.PPI", "STRING")
  
  ks_files = purrr::map(ks_libs, function(in_pattern){
    grep(in_pattern, all_files, value = TRUE)
  })
  ppi_files = purrr::map(ppi_libs, function(in_pattern){
    grep(in_pattern, all_files, value = TRUE)
  })
  
  ks_len = purrr::map_int(ks_files, length)
  ppi_len = purrr::map_int(ppi_files, length)
  ks_files = unlist(ks_files[ks_len > 0])
  ppi_files = unlist(ppi_files[ppi_len > 0])
  
  ks_data = purrr::map(ks_files, parse_gmt)
  ppi_data = purrr::map(ppi_files, parse_gmt)
  
  ks_names = purrr::map(ks_data, ~ names(.x))
  all_ks = unique(c(unlist(ks_data), unlist(ks_names)))
  ppi_names = purrr::map(ppi_data, ~ names(.x))
  all_ppi = unique(c(unlist(ppi_data), unlist(ppi_names)))
  
  
  ks_matrix = matrix(0, nrow = length(all_ks), ncol = length(all_ks))
  rownames(ks_matrix) = colnames(ks_matrix) = all_ks
  for (iks in ks_data) {
    for (ks_int in names(iks)) {
      #message(ks_int)
      ks_matrix[iks[[ks_int]], ks_int] = ks_matrix[iks[[ks_int]], ks_int] + 1
    }
  }
  
  ppi_matrix = matrix(0, nrow = length(all_ppi), ncol = length(all_ppi))
  rownames(ppi_matrix) = colnames(ppi_matrix) = all_ppi
  for (ippi in ppi_data) {
    for (ppi_int in names(ippi)) {
      #message(ppi_int)
      if (ppi_int %in% rownames(ks_matrix)){
        ks_tmp = ks_matrix[ppi_int, ]
        ks_match = intersect(names(ks_tmp)[ks_tmp > 0], ippi[[ppi_int]])
        ppi_novel = setdiff(ippi[[ppi_int]], names(ks_tmp)[ks_tmp > 0])
        ks_matrix[ks_match, ppi_int] = ks_matrix[ks_match, ppi_int] + 1
        
      } else {
        ppi_novel = ippi[[ppi_int]]
      }
      ppi_matrix[ppi_novel, ppi_int] = ppi_matrix[ppi_novel, ppi_int] + 1
    }
  }
  
  
  ks_df = purrr::map_dfr(colnames(ks_matrix), function(col_name){
    tmp_col = ks_matrix[, col_name]
    tmp_col = tmp_col[tmp_col > 0]
    if (length(tmp_col) > 0) {
      out_df = data.frame(from = col_name, to = names(tmp_col), weight = tmp_col)
      return(out_df)
    } else {
      return(NULL)
    }
  })
  ks_df$type = "kinase-substrate"
  
  ppi_df = purrr::map_dfr(colnames(ppi_matrix), function(col_name){
    tmp_col = ppi_matrix[, col_name]
    tmp_col = tmp_col[tmp_col > 0]
    if (length(tmp_col) > 0) {
      out_df = data.frame(from = col_name, to = names(tmp_col), weight = tmp_col)
      return(out_df)
    } else {
      return(NULL)
    }
  })
  ppi_df2 = data.frame(from = ppi_df$to,
                       to = ppi_df$from,
                       weight = ppi_df$weight)
  ppi_df$type = "ppi"
  ppi_df2$type = "ppi"
  ppi_full = rbind(ppi_df, ppi_df2)
  ks_graph = tidygraph::as_tbl_graph(ks_df)
  ppi_graph = tidygraph::as_tbl_graph(ppi_df)
  full_graph = tidygraph::graph_join(ks_graph, ppi_graph)
  
  all_kinases = unique(c(unlist(ks_names), unlist(ppi_names)))
  full_graph = full_graph %>%
    activate(nodes) %>%
    mutate(type = case_when(
      name %in% all_kinases ~ "kinase",
      TRUE ~ "other"
    ))
  return(full_graph)
}

get_lfc = function(lfc_files){
  all_lfc = purrr::map_dfr(lfc_files, function(in_file){
    tmp_df = read.table(in_file, header = TRUE, sep = "\t",
                        stringsAsFactors = FALSE)
    if (grepl("ptk", in_file)) {
      source = "PTK"
    } else {
      source = "STK"
    }
    tmp_df2 = tmp_df %>%
      dplyr::mutate(peptide2 = gsub("_.*", "", Peptide)) %>%
      dplyr::group_by(peptide2) %>%
      dplyr::summarise(protein_mean_lfc = mean(LFC))
    tmp_df2$source = source
    tmp_df2
  })
  
  dup_prots = all_lfc %>%
    dplyr::group_by(peptide2) %>%
    dplyr::summarise(n_prot = n()) %>%
    dplyr::filter(n_prot > 1)
  
  bad_dups = all_lfc %>%
    dplyr::filter(peptide2 %in% dup_prots$peptide2) %>%
    dplyr::group_by(peptide2) %>%
    dplyr::summarise(lfc_ratio = protein_mean_lfc[1] / protein_mean_lfc[2]) %>%
    dplyr::filter(lfc_ratio < 0)
  
  good_dups = all_lfc %>%
    dplyr::filter(peptide2 %in% dup_prots$peptide2) %>%
    dplyr::group_by(peptide2) %>%
    dplyr::summarise(lfc_ratio = protein_mean_lfc[1] / protein_mean_lfc[2]) %>%
    dplyr::filter(lfc_ratio > 0)
  
  dup_2_single = all_lfc %>%
    dplyr::filter(peptide2 %in% good_dups$peptide2) %>%
    dplyr::group_by(peptide2) %>%
    dplyr::slice(n = 1)
  
  all_lfc = all_lfc %>%
    dplyr::filter(!(peptide2 %in% dup_prots$peptide2))
  
  all_lfc = rbind(all_lfc, dup_2_single)
  
  all_lfc = all_lfc %>%
    dplyr::mutate(direction = dplyr::case_when(
      protein_mean_lfc < 0 ~ "neg",
      protein_mean_lfc > 0 ~ "pos"
    ))
  
  all_lfc
    
}

create_enriched_network = function(kea3_results, n_gene = 10, kea3_network, tissue_list){
  genes_enrich = purrr::map_dfr(kea3_results, function(.x){
    data.frame(gene = .x$`Integrated--meanRank`$TF[seq_len(n_gene)],
               type = "kinase",
               enrichment = "enriched")
  })
  genes_enrich = genes_enrich %>%
    dplyr::filter(gene %in% tissue_list)
  
  # and overlapping genes, essentially the kinase targets
  genes_overlaps = purrr::map_dfr(kea3_results,
                                  function(.x){
                                    overlap_genes = .x$`Integrated--meanRank`$Overlapping_Genes[seq_len(n_gene)]
                                    all_over = unlist(strsplit(overlap_genes, ","))
                                    data.frame(gene = all_over,
                                               type = "target",
                                               enrichment = "overlap")
                                  }) %>%
    unique()
  rownames(genes_overlaps) = NULL
  genes_overlaps = genes_overlaps %>%
    dplyr::filter(gene %in% tissue_list)
  
  genes_overlaps = genes_overlaps %>%
    dplyr::filter(!(gene %in% genes_enrich$gene))
  
  enrich_results = rbind(genes_enrich, genes_overlaps) %>%
    dplyr::filter(gene %in% tissue_list)
  
  enrich_network = kea3_network %>%
    activate(nodes) %>%
    filter(name %in% enrich_results$gene) %>%
    mutate(enriched = case_when(
      name %in% genes_enrich$gene ~ "enriched",
      TRUE ~ "overlap"
    )) %>%
    activate(edges) %>%
    filter(weight > 1) %>%
    activate(nodes) %>%
    mutate(degree = local_size(order = 1, mindist = 1)) %>% # how many neighbors?
    filter(degree > 1)
  
  enrich_network_ks = enrich_network %>%
    activate(edges) %>%
    filter(type %in% "kinase-substrate") %>%
    activate(nodes) %>%
    mutate(degree = local_size(order = 1, mindist = 1)) %>%
    filter(degree > 1) %>%
    mutate(type2 = case_when(
      (type %in% "kinase") & (enriched %in% "enriched") ~ "kinase.E",
      (type %in% "kinase") & (enriched %in% "overlap") ~ "kinase.O",
      TRUE ~ "other"
    ))
  
  node_info = enrich_network_ks %>%
    activate(nodes) %>%
    data.frame() %>%
    dplyr::rename(TypeSource = type2)
  
  list(genes = node_info, network = enrich_network_ks)
  
}

remove_nodes_other = function(full_network, keep_genes){
  
}

annotate_full_network = function(full_network, network_1, network_2,
                                 network_1_id = "human", network_2_id = "mouse"){
  
  full_nodes = full_network %>%
    activate(nodes) %>%
    data.frame()
  full_edges = full_network %>%
    activate(edges) %>%
    data.frame()
  full_edges$from = full_nodes$name[full_edges$from]
  full_edges$to = full_nodes$name[full_edges$to]
  full_edges$from_to = paste0(full_edges$from, ".", full_edges$to)
  
  
  n1_edges = network_1 %>%
    activate(edges) %>%
    data.frame()
  n1_nodes = network_1 %>%
    activate(nodes) %>%
    data.frame()
  n1_edges$from = n1_nodes$name[n1_edges$from]
  n1_edges$to = n1_nodes$name[n1_edges$to]
  n1_edges$from_to = paste0(n1_edges$from, ".", n1_edges$to)
  
  n2_edges = network_2 %>%
    activate(edges) %>%
    data.frame()
  n2_nodes = network_2 %>%
    activate(nodes) %>%
    data.frame()
  n2_edges$from = n2_nodes$name[n2_edges$from]
  n2_edges$to = n2_nodes$name[n2_edges$to]
  n2_edges$from_to = paste0(n2_edges$from, ".", n2_edges$to)
  
  full_edges = full_edges %>%
    dplyr::mutate(network = dplyr::case_when(
      (from_to %in% n1_edges$from_to) & (from_to %in% n2_edges$from_to) ~ "both",
      from_to %in% n1_edges$from_to ~ network_1_id,
      from_to %in% n2_edges$from_to ~ network_2_id,
    ))
  
  new_network = tidygraph::as_tbl_graph(full_edges)
  new_network = new_network %>%
    activate(nodes) %>%
    mutate(network = case_when(
      (name %in% n1_nodes$name) & (name %in% n2_nodes$name) ~ "both",
      name %in% n1_nodes$name ~ network_1_id,
      name %in% n2_nodes$name ~ network_2_id
    ))
  
  new_network
  
}

subset_kegg_network = function(full_network, subset_nodes){
  node_df = full_network %>%
    activate(nodes) %>%
    data.frame()
  edge_df = full_network %>%
    activate(edges) %>%
    data.frame()
  edge_named_df = data.frame(from = node_df$name[edge_df$from],
                             to = node_df$name[edge_df$to])
}

subset_network_alpha = function(full_network, subset_results, keep_network = c("both", "human"), exclude_alpha = 0){
  # for testing the function and iterating through it
  #full_network = all_results_source
  #subset_results = human_results
  #keep_network = c("both", "human")
  #exclude_alpha = 0
  
  subset_source = split(subset_results$genes$name, subset_results$genes$TypeSource)
  
  subset_network = full_network %>%
    activate(edges) %>%
    mutate(alpha = case_when(
      network %in% keep_network ~ 1,
      TRUE ~ exclude_alpha
    )) %>%
    activate(nodes) %>%
    mutate(alpha = case_when(
      network %in% keep_network ~ 1,
      TRUE ~ exclude_alpha
    ),
    org_name = name,
    name = case_when(
      alpha == 1 ~ org_name,
      TRUE ~ ""
    ))
  
  subset_network = subset_network %>%
    activate(nodes) %>%
    mutate(TypeSource = as.factor(case_when(
      name %in% subset_source$kinase.E ~ "kinase.E",
      name %in% subset_source$kinase.O ~ "kinase.O",
      name %in% subset_source$other ~ "other",
      TRUE ~ "NA"
    )))
  subset_network
}

export_plots = function(plot_list, ppt_file){
  new_ppt = officer::read_pptx()
  for (iplot in plot_list) {
    new_ppt = officer::add_slide(new_ppt, layout = "Blank")
    if (inherits(iplot, "ggplot")) {
      p_dml = rvg::dml(ggobj = iplot)
      
      new_ppt = officer::ph_with(new_ppt, value = p_dml, location = officer::ph_location(width = 8, height = 8))
      
    } else if (inherits(iplot, "Heatmap")) {
      
      p_dml = rvg::dml(code = draw(iplot, merge_legend = TRUE))
      
      new_ppt = officer::ph_with(new_ppt, value = p_dml, location = officer::ph_location(width = 8, height = 8))
    
    } else if (inherits(iplot, "flextable")) {
      new_ppt = officer::ph_with(new_ppt, value = iplot, location = officer::ph_location_left())
    }
    
  }
  print(new_ppt, target = here::here("data", "outputs", "figures", ppt_file))
}

get_kinase_lfc = function(){
  mapping = c("Q" = "\\[theta\\]", "D" = "\\[delta\\]", "G" = "\\[gamma\\]", 
              "H" = "\\[eta\\]", "I" = "\\[iota\\]", "E" = "\\[epsilon\\]", 
              "Z" = "\\[zeta\\]", "A" = "\\[alpha\\]", "B" = "\\[beta\\]")
  lfc_dir = here::here("data", "inputs", "2021_10_07-Justin-Data-Deposit", "UKA derived output fold-changes")
  lfc_files = dir(lfc_dir, pattern = "txt$", full.names = TRUE)
  
  lfc_values = purrr::map_df(lfc_files, function(.x){
    df = read.table(.x, header = TRUE, sep = "\t")
    if (grepl("human", .x)) {
      df$organism = "Human"
    } else {
      df$organism = "Mouse"
    }
    
    if (grepl("PTK", .x)) {
      df$type = "PTK"
    } else {
      df$type = "STK"
    }
    df
  })
  lfc_values$KN2 = lfc_values$Kinase.Name
  for (imapping in names(mapping)) {
    has_map = grepl(mapping[imapping], lfc_values$Kinase.Name)
    lfc_values[has_map, "KN2"] = gsub(mapping[imapping], imapping, lfc_values[has_map, "Kinase.Name"])
  }
  lfc_values$name = toupper(lfc_values$KN2)
  lfc_values$LFC = lfc_values$Mean.Final.Score
  split_lfc = split(lfc_values[, c("name", "LFC")], lfc_values$organism)
  lfc_out = dplyr::left_join(split_lfc$Mouse, split_lfc$Human, by = "name",
                             suffix = c(".Mouse", ".Human"))
  lfc_out
}

# given a graph, a set of nodes, and a minimum edge weight
# to filter out low confidence edges, subset to the provided
# nodes and their immediate neighbors.
subset_nodes_neighbors = function(in_graph, nodes = NULL, edge_min = 1, edge_type = "kinase-substrate"){
  if (is.null(nodes)) {
    stop("nodes = NULL\nYou didn't provide which nodes you want to start with!")
  }
  
  filter_edges = in_graph %>%
    activate(edges) %>%
    filter(weight >= edge_min, type %in% edge_type) %>%
    remove_stranded_nodes()
  
  all_nodes = filter_edges %>%
    activate(nodes) %>%
    data.frame()
  
  node_loc = which(all_nodes$name %in% nodes)
  
  neighbor_edges = filter_edges %>%
    activate(edges) %>%
    filter((from %in% node_loc) | (to %in% node_loc)) %>%
    remove_stranded_nodes()
  
  neighbor_edges
}

remove_stranded_nodes = function(in_graph){
  in_graph %>%
    activate(nodes) %>%
    mutate(degree = local_size(order = 1, mindist = 1)) %>%
    filter(degree > 0)
}

read_raw_files = function(){
  near_raw_dir = here::here("data", "inputs", "pamgene_near_raw")
  raw_files = dir(near_raw_dir, pattern = ".txt", full.names = TRUE)
  
  in_raw = function(raw_file){
    raw_data = scan(raw_file, what = character(), sep = "\n", quiet = TRUE)
    raw_data = raw_data[-1]
    id_loc = raw_data %in% "ID"
    raw_data = raw_data[!id_loc]
    all_data = purrr::map_dfc(raw_data, function(raw_chr){
      tmp_chr = strsplit(raw_chr, "\t")[[1]]
      
      top_chr = tmp_chr[c(1, 2)]
      tmp_chr = tmp_chr[c(-1, -2)]
      var_name = top_chr[nchar(top_chr) > 0]
      
      #message(var_name)
      
      if (!is.na(suppressWarnings(as.numeric(tmp_chr[1])))) {
        out_var = as.numeric(tmp_chr)
      } else {
        out_var = tmp_chr
      }
      out_df = data.frame(v = out_var)
      names(out_df) = var_name
      out_df
    })
    all_data
  }
  
  raw_data = purrr::map(raw_files, in_raw)
  names(raw_data) = basename(raw_files)
  raw_data
}

transform_raw_data = function(raw_data){
  raw_full_info = data.frame(barcode = raw_data$Barcode,
                             array = raw_data$Array,
                             exposure = raw_data$`Exposure time`,
                             cycle = raw_data$Cycle,
                             sample_name = raw_data$`Sample name`,
                             comment = raw_data$`Comment 1`) %>%
    dplyr::mutate(replicate = case_when(
      grepl("3$", sample_name) ~ "3",
      grepl("2$", sample_name) ~ "2",
      grepl("1$", sample_name) ~ "1"
    ))
  
  raw_info = raw_full_info %>%
    dplyr::select(barcode, array, sample_name, comment, replicate) %>%
    unique()
  last_ref = max(which(grepl("REF.*", names(raw_data))))
  start_data = last_ref + 1
  end_data = ncol(raw_data)
  raw_vals = raw_data[, seq(start_data, end_data)]
  raw_info_vals = cbind(raw_full_info, raw_vals)
  
  raw_long_values = raw_info_vals %>%
    tidyr::pivot_longer(c(-barcode, -array, -exposure, -cycle, -sample_name, -comment, -replicate), names_to = "peptide")
  raw_long_values = raw_long_values %>%
    dplyr::mutate(measure_id = paste0(peptide, ".", exposure, ".", cycle), sample_name2 = paste0(comment, "_", replicate))
  list(data = raw_long_values, sample_info = raw_info)
}

create_go_annotation = function(db, ontology = NULL){
  all_genes = keys(db)
  go_all_gene = AnnotationDbi::select(db, keys = all_genes, columns = c("GOALL", "ONTOLOGYALL"))
  
  if (!is.null(ontology)) {
    go_all_gene = go_all_gene[go_all_gene$ONTOLOGYALL == ontology, ]
    ontology_type = paste0("GO.", ontology)
  } else {
    ontology_type = "GO.all"
  }
  go_2_gene = split(go_all_gene$ENTREZID, go_all_gene$GOALL)
  go_2_gene = lapply(go_2_gene, unique)
  go_desc = AnnotationDbi::select(GO.db::GO.db, keys = names(go_2_gene), columns = "TERM", keytype = "GOID")$TERM
  names(go_desc) = names(go_2_gene)

  go_annotation = categoryCompare2::annotation(annotation_features = go_2_gene,
                                               description = go_desc,
                                               annotation_type = ontology_type,
                                               feature_type = "ENTREZID")
  go_annotation
}


create_kegg_annotation = function(kegg_data){
  kegg_data$description = gsub(" - Homo sapiens (human)", "", kegg_data$description, fixed = TRUE)
  
  kegg_annotation = categoryCompare2::annotation(
    annotation_features = kegg_data$annotation,
    description = kegg_data$description[names(kegg_data$annotation)],
    annotation_type = "kegg",
    feature_type = "ENTREZID"
  )
  kegg_annotation
}


extract_stats_table = function(in_results, specific_entrez){
  sig_res = in_results@statistics@significant@significant
  which_sig = apply(sig_res, 1, function(.x){
    sum(.x > 0)
  }) > 0
  stat_table = in_results@statistics@statistic_data[which_sig, ]
  stat_table = stat_table %>%
    tibble::rownames_to_column("path_id") %>%
    mutate(pathway = in_results@annotation@description[path_id])
  stat_table = add_specific_count(stat_table,
                                  in_results@annotation,
                                  specific_entrez)
  stat_table$n_annotated = in_results@annotation@counts[stat_table$path_id]
  stat_table
}

add_specific_count = function(results_table, annotation_obj, specific_entrez){
  n_spec = purrr::imap_dfr(annotation_obj@annotation_features,
                          function(.x, .y){
    data.frame(path_id = .y, n_specific = sum(specific_entrez %in% .x))
                          })
  results_table = dplyr::left_join(results_table, n_spec, by = "path_id")
  results_table
}

add_genes = function(results_table, sig_obj, gene_df){
  use_ids = results_table$path_id
  all_features = union(sig_obj@enriched[[1]]@features, sig_obj@enriched[[2]]@features)
  annot_obj = sig_obj@annotation@annotation_features
  annotation_genes = purrr::map_chr(use_ids, function(in_id){
    tmp_features = intersect(all_features, annot_obj[[in_id]])
    gene_df %>%
      dplyr::filter(ENTREZID %in% tmp_features) %>%
      dplyr::pull(SYMBOL) %>%
      sort() %>%
      paste(., collapse = ":")
  })
  results_table$genes = annotation_genes
  results_table
}

export_pathway_tables = function(table_list, pathway_file){
  openxlsx::write.xlsx(table_list, file = pathway_file, overwrite = TRUE)
}

score_annotation_overlap = function(annot_list){
  overlap_matrix = matrix(0, nrow = length(annot_list), ncol = length(annot_list))
  colnames(overlap_matrix) = rownames(overlap_matrix) = names(annot_list)
  
  for (i in seq(1, length(annot_list))) {
    for (j in seq(i, length(annot_list))) {
      overlap_matrix[i, j] = categoryCompare2::combined_coefficient(annot_list[[i]], annot_list[[j]])
    }
  }
  dist_overlap = as.dist(1 - overlap_matrix)
  annot_hclust = as.dendrogram(hclust(dist_overlap))
  #order_annot = dendsort::dendsort(annot_hclust, type = "average")
  annot_hclust
}