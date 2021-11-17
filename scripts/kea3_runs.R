# this is the code for querying the kea3 api.
# We have it in an R script so that we have to choose to run it,
# and don't have to requery the API each time we want to update a report.

kea3_outfile = here::here("data", "outputs", "rds_files", "kea3_results.rds")

# get the list of files
kea3_input_lists = dir(here::here("data", "inputs", "kea3"), pattern = "^z", full.names = TRUE)

# for each file, read in the list of genes and remove the header
read_genes = function(gene_file){
  scan(gene_file, character(), skip = 1, sep = "\n")
}

# apply the read_genes function to each file
kea3_input_genes = purrr::map(kea3_input_lists, read_genes)
names(kea3_input_genes) = basename(kea3_input_lists)

# query the api
query_kea3 = function(gene_list, query_name){
  url = "https://maayanlab.cloud/kea3/api/enrich/"
  input_payload = list(query_name = query_name, gene_set = gene_list)
  
  api_response = httr::POST(url = url, body = input_payload, encode = "json")
  json_output = httr::content(api_response, "text")
  jsonlite::fromJSON(json_output)
}

# we use "imap" here to pass the list of genes,
# and the name of the entry
kea3_results = purrr::imap(kea3_input_genes, query_kea3)

saveRDS(kea3_results, kea3_outfile)