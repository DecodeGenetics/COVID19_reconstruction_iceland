library(igraph)
library(ggraph)
library(readr)
library(tidyr)
library(dplyr)
library(readr)

path_results <- ''
model_name <- ''
outbreak_dat <- read_tsv(paste0(path_results,'outbreak_dat/', model_name, '.tsv'))
outbreaker_mod <- readRDS(paste0(path_results,'outbreaker_res/', model_name, '.rds'))

N <- 2522

tree_no <- 199
alpha <- c()
idx <- 1:N
for (i in idx) {
  alpha[i] <- outbreaker_mod[[paste0('alpha_', i)]][tree_no]
}
root <- which(is.na(alpha))

idx <- idx[!is.na(alpha)]
alpha <- alpha[!is.na(alpha)]
leaves <- idx[!idx %in% alpha]

in_q <- outbreak_dat %>% filter(q_type_binary == 'in_q') %>% .$idx

vc <- rep('Outside quarantine',N)
vc[in_q] <- 'In quarantine'
vc[root] <- 'Root'

el <- matrix(c(alpha, idx), nc=2)
G <- graph_from_edgelist(el)
tiff(filename='infection_tree.tiff', width=1870, height=2048, units='px',bg='#C6E8DA')
ggraph(G, layout='tree', circular=T) +
  geom_edge_link(arrow = arrow(length = unit(5, 'mm')), end_cap = circle(5, 'mm'), edge_width=1.5) +
  geom_node_point(size = 10, shape=21, color='black', stroke=2, aes(fill = as.factor(vc))) +
  scale_fill_manual(values=c('#03adfc', '#fc0303', '#fcba03')) +
  theme_void() +
  theme(legend.position='none')
dev.off()
                    
