#### ---- Packages ---- ####
library(dplyr)
library(data.table)
library(ggplot2)
library(gridExtra)
library(tidyr)
library(scales)
library(readr)
library(here)
library(ggsci)
source(here('help_functions_analysis.R'))
source(here('help_functions_immunity.R'))
source(here('help_functions_immunity_analysis.R'))

###############################################   Vaccination figure for paper #######################################
out_dir <- 'results/immunity_simulation/'
fig_out_dir <- 'figures/immunty_simulation/'

immunity_res_prior <- read_immunity_files(subtrees=F,replay=T,init_size=0,efficacy=0.6,N=50,max_chains=1,max_sim=19,path=out_dir)
immunity_res_prior_90p <- read_immunity_files(subtrees=F,replay=T,init_size=0,efficacy=0.9,N=50,max_chains=1,max_sim=19,path=out_dir)
immunity_res <- read_immunity_files(subtrees=F,replay=F,init_size=0,efficacy=0.6,N=50,max_chains=1,max_sim=19,path=out_dir)
immunity_res_90p <- read_immunity_files(subtrees=F,replay=F,init_size=0,efficacy=0.9,N=50,max_chains=1,max_sim=19,path=out_dir)

write.table(immunity_res_prior,file=paste0(out_dir,'/immunity_res_prior.tsv'),row.names=F,quote=F,sep='\t')
write.table(immunity_res_prior_90p,file=paste0(out_dir,'/immunity_res_prior_90p.tsv'),row.names=F,quote=F,sep='\t')
write.table(immunity_res,file=paste0(out_dir,'/immunity_res.tsv'),row.names=F,quote=F,sep='\t')
write.table(immunity_res_90p,file=paste0(out_dir,'/immunity_res_90p.tsv'),row.names=F,quote=F,sep='\t')

immunity_res_50_prior <- read_immunity_files(subtrees=F,replay=T,init_size=50,efficacy=0.6,N=50,max_chains=1,max_sim=19,path=out_dir)
immunity_res_50_prior_90p <- read_immunity_files(subtrees=F,replay=T,init_size=50,efficacy=0.9,N=50,max_chains=1,max_sim=19,path=out_dir)
immunity_res_50 <- read_immunity_files(subtrees=F,replay=F,init_size=50,efficacy=0.6,N=50,max_chains=1,max_sim=19,path=out_dir)
immunity_res_50_90p <- read_immunity_files(subtrees=F,replay=F,init_size=50,efficacy=0.9,N=50,max_chains=1,max_sim=19,path=out_dir)

immunity_res_subtrees_prior <- read_immunity_files(subtrees=T,replay=T,init_size=0,efficacy=0.6,N=50,max_chains=1,max_sim=19,path=out_dir)
immunity_res_subtrees_prior_90p <- read_immunity_files(subtrees=T,replay=T,init_size=0,efficacy=0.9,N=50,max_chains=1,max_sim=19,path=out_dir)
immunity_res_subtrees <- read_immunity_files(subtrees=T,replay=F,init_size=0,efficacy=0.6,N=50,max_chains=1,max_sim=19,path=out_dir)
immunity_res_subtrees_90p <- read_immunity_files(subtrees=T,replay=F,init_size=0,efficacy=0.9,N=50,max_chains=1,max_sim=19,path=out_dir)

figure3 <- grid.arrange(plot_immunity_prior(immunity_res_prior) + xlab('') + theme(legend.position = 'none') + ggtitle('A. Prior vaccinations - First dose'),
                        plot_immunity(immunity_res) + xlab('') + ylab('') + ggtitle('B. First dose'),
                        plot_immunity_prior(immunity_res_prior_90p) + theme(legend.position = 'none')+ ggtitle('C. Prior vaccinations - Second dose'),
                        plot_immunity(immunity_res_90p) + ylab('') + theme(legend.position = 'none') + ggtitle('D. Second dose'),ncol=2)
ggsave(figure3,filename = paste0(fig_out_dir,'/figure3.png'),width=180,height=180,units = 'mm',dpi = 320)

sup_figure_50 <- grid.arrange(plot_immunity_prior(immunity_res_50_prior) + xlab('') + theme(legend.position = 'none') + ggtitle('A. Prior vaccinations - First dose'),
                              plot_immunity(immunity_res_50) + xlab('') + ylab('') + ggtitle('B. First dose'),
                              plot_immunity_prior(immunity_res_50_prior_90p) + theme(legend.position = 'none')+ ggtitle('C. Prior vaccinations - Second dose'),
                              plot_immunity(immunity_res_50_90p) + ylab('') + theme(legend.position = 'none') + ggtitle('D. Second dose'),ncol=2)
ggsave(sup_figure_50,filename = paste0(fig_out_dir,'/sup_figure_50.png'),width=180,height=180,units = 'mm',dpi = 320)

sup_figure_prop <- grid.arrange(plot_immunity_subtrees_prior(immunity_res_subtrees_prior) + xlab('') + theme(legend.position = 'none') + ggtitle('A. Prior vaccinations - First dose'),
                                plot_immunity_subtrees(immunity_res_subtrees) + xlab('') + ylab('') + ggtitle('B. First dose'),
                                plot_immunity_subtrees_prior(immunity_res_subtrees_prior_90p) + theme(legend.position = 'none')+ ggtitle('C. Prior vaccinations - Second dose'),
                                plot_immunity_subtrees(immunity_res_subtrees_90p) + ylab('') + theme(legend.position = 'none') + ggtitle('D. Second dose'),ncol=2)
ggsave(sup_figure_prop,filename = paste0(fig_out_dir,'/sup_figure_prop.png'),width=180,height=180,units = 'mm',dpi = 320)

##### Visualize subtree roots
path_results <- 'results/full_model/'
model_names <- c("chain_1","chain_2",
                 "chain_3","chain_4")
burnin <- 5000
outbreak_dat <- read_tsv(paste0(path_results,'outbreak_dat/',model_names[1],'.tsv')) %>%
  mutate(shows_symptoms=!is.na(date_symptoms))
outbreaker_obj <- readRDS(paste0(path_results,'outbreaker_obj/',model_names[1],'.rds'))

ctd_dat <- as_tibble(outbreaker_obj$ctd) %>%
  dplyr::select(idx_from=idx_to,idx_to=idx_from) %>%
  bind_rows(as_tibble(outbreaker_obj$ctd)) %>%
  inner_join(dplyr::select(outbreak_dat,t_id,idx),by=c('idx_to'='idx')) %>%
  inner_join(dplyr::select(outbreak_dat,t_id,idx),by=c('idx_from'='idx'),suffix=c('_to','_from')) %>%
  dplyr::select(t_id_from,t_id_to) %>%
  mutate(contact=T) %>%
  distinct()

num_trees_per_chain <- (10000-burnin)/50 + 1
outbreaker_res_dat <- lapply(1:4,function(i){
  outbreaker_mod <- readRDS(paste0(path_results,'outbreaker_res/',model_names[i],'.rds'))
  get_outbreaker_res_long(outbreaker_mod,outbreak_dat,burnin=burnin) %>% 
    inner_join(dplyr::select(outbreak_dat,t_id,date,date_q,date_exposure,date_symptoms,household_id),by='t_id') %>%
    left_join(select(outbreak_dat,t_id_from=t_id,household_id),by='t_id_from',suffix=c('','_from')) %>%
    left_join(ctd_dat,by=c('t_id'='t_id_to','t_id_from')) %>%
    mutate(same_household=if_else(is.na(household_id) | is.na(household_id_from),F,household_id==household_id_from),
           contact=if_else(is.na(contact),F,contact),
           chain=i,
           iter=(i-1)*num_trees_per_chain+iter) %>%
    select(chain,everything())
}) %>% bind_rows()

outdegree_dat <- get_outbreaker_outdegree_dat(outbreaker_res_dat,outbreak_dat) %>%
                  inner_join(dplyr::select(outbreak_dat,t_id,ct,date,shows_symptoms,score,age_at_diag,child,q_type_binary,q_type),by='t_id') %>%
                  unite(col = child_q,child,q_type_binary,sep = ',',remove = F)
trees <- create_trees(outbreaker_res_dat)
subtrees <- get_subtree_size(outbreaker_res_dat,trees,50) %>%
            inner_join(select(outbreak_dat,t_id,idx,age_at_diag),by=c('root'='idx'))

ggplot(filter(subtrees,n>=100),aes(age_at_diag)) +
geom_histogram(bins=20,col='black',fill=my_colors(2)[2],size=0.3) +
theme_classic() + 
immunity_theme() +
  theme(axis.title.y=element_text(size=8),
        axis.title.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.text.x=element_text(size=8),
        axis.line = element_line(colour = 'black', size = 0.3))+
ylab('Number of subtree roots') +
xlab('Age')
ggsave(filename = paste0(fig_out_dir,'/root_age_distr.png'),width=90,height=60,units = 'mm',dpi = 320)

ggplot(subtrees,aes(n)) +
geom_histogram(bins=48,col='black',fill=my_colors(2)[2],size=0.3) +
geom_vline(xintercept=100,col='red',linetype='dashed',size=0.3) +
theme_classic() + 
immunity_theme() +
theme(axis.title.y=element_text(size=8),
      axis.title.x=element_text(size=8),
      axis.text.y=element_text(size=8),
      axis.text.x=element_text(size=8),
      axis.line = element_line(colour = 'black', size = 0.3)) +
scale_y_continuous(trans=scales::log10_trans(),breaks = trans_breaks("log10", function(x) 10^x)) +
ylab('Number of subtrees') +
xlab('Subtree size')
ggsave(filename = paste0(fig_out_dir,'/subtree_size_distr.png'),width=90,height=60,units = 'mm',dpi = 320)

#distribution of outdegrees of roots
inner_join(subtrees,select(outdegree_dat,iter,t_id,outdegree_all),by=c('iter','t_id')) %>%
  mutate(size_cat=cut(n,breaks=c(-Inf,100,Inf),labels=c('Subtree of size < 100','Subtree of size >= 100'))) %>%
ggplot(aes(size_cat,outdegree_all)) +
  geom_boxplot() +
  theme_classic() + 
  immunity_theme() +
  theme(axis.title.y=element_text(size=8),
        axis.title.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.text.x=element_text(size=8),
        axis.line = element_line(colour = 'black', size = 0.3)) +
  ylab('Out-degree') +
  xlab('Subtree size')
ggsave(filename = paste0(fig_out_dir,'/subtree_size_distr.png'),width=90,height=60,units = 'mm',dpi = 320)

