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
  
  ppi_matrix = matrix(0, nrow = length(all_ppi), ncol = length(all_ppi))
  rownames(ppi_matrix) = colnames(ppi_matrix) = all_ppi
  for (ippi in ppi_data) {
    for (ppi_int in names(ippi)) {
      #message(ppi_int)
      ppi_matrix[ppi_int, ippi[[ppi_int]]] = ppi_matrix[ppi_int, ippi[[ppi_int]]] + 1
    }
  }
  ks_graph = tidygraph::as_tbl_graph(ks_matrix, mode = "out")
  ppi_graph = tidygraph::as_tbl_graph(ppi_matrix)
}

parse_gmt = function(gmt_file){
  gmt = scan(gmt_file, what = "", sep = "\n", quiet = TRUE)
  gmt_split = strsplit(gmt, "\t")
  names(gmt_split) = purrr::map_chr(gmt_split, ~ gsub("_.*", "", .x[[1]]))
  gmt_split = purrr::map(gmt_split, ~ .x[-1])
  gmt_split
}

create_interaction = function()