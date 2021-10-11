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
  kea3_data_dir = here::here("kea3_datasets")
  all_files = dir(kea3_data_dir, full.names = TRUE)
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
  
  ppi_matrix = matrix(0, nrow = length(all_ppi), ncol = length(all_ppi))
  rownames(ppi_matrix) = colnames(ppi_matrix) = all_ppi
  for (ippi in ppi_data) {
    for (ppi_int in names(ippi)) {
      #message(ppi_int)
      ppi_matrix[ppi_int, ippi[[ppi_int]]] = ppi_matrix[ppi_int, ippi[[ppi_int]]] + 1
    }
  }
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
    p_dml = rvg::dml(ggobj = iplot)
    new_ppt = officer::add_slide(new_ppt, layout = "Blank")
    new_ppt = officer::ph_with(new_ppt, value = p_dml, location = officer::ph_location(width = 8, height = 8))
  }
  print(new_ppt, target = here::here(ppt_file))
}

get_kinase_lfc = function(){
  mapping = c("Q" = "\\[theta\\]", "D" = "\\[delta\\]", "G" = "\\[gamma\\]", 
              "H" = "\\[eta\\]", "I" = "\\[iota\\]", "E" = "\\[epsilon\\]", 
              "Z" = "\\[zeta\\]", "A" = "\\[alpha\\]", "B" = "\\[beta\\]")
  lfc_dir = here::here("input/2021_10_07-Justin-Data-Deposit/UKA derived output fold-changes")
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