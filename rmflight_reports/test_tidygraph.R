library(tidygraph)
library(ggraph)
tmp_graph = create_ring(10)
tmp_graph = tmp_graph %>%
  activate(nodes) %>%
  mutate(name = seq(1, 10))

ggraph(tmp_graph, layout = "graphopt") +
  geom_edge_link() +
  geom_node_label(aes(label = name))

tmp_graph2 = tmp_graph %>%
  activate(nodes) %>%
  filter(!(name %in% 5))
ggraph(tmp_graph2, layout = "graphopt") +
  geom_edge_link() +
  geom_node_label(aes(label = name))
