#### ---- Packages ---- ####
library(argparser)
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(tidyr)
library(readxl)
library(here)
source(here('help_functions_analysis.R'))
source(here('help_functions_immunity.R'))

#### Input arguments ####
p <- arg_parser("Simulate vaccination strategies","immunity_infection_trees.R")
p <- add_argument(p,"--model_name",nargs=1,type='character',
                  help="Name of outbreaker model")
p <- add_argument(p,"--sim_nr",nargs=1,type='numeric',
                  help="Incremental number of simulation")
p <- add_argument(p,"-N",nargs=1,type='numeric',default=50,
                  help="Number of simulations to perform.")
p <- add_argument(p,"--first_dose_efficacy",nargs=1,type='numeric',default=0.6,
                  help="Desired efficacy of first dose")
p <- add_argument(p,"--second_dose_efficacy",nargs=1,type='numeric',default=0.9,
                  help="Desired efficacy of second dose")
p <- add_argument(p,"--init_size",nargs=1,type='numeric',default=0,
                  help="Size of initial outbreak")
p <- add_argument(p,"--step_size",nargs=1,type='numeric',default=0.05,
                  help="Percentage interval at which to simulate")
p <- add_argument(p,"--replay",flag=T,
                  help="Indicates whether we should replay prior vaccinations or not")
p <- add_argument(p,"--subtrees",flag=T,
                  help="Indicates whether we should simulate subtrees or not")
p <- add_argument(p,"--out_dir",type='character',
                  help="path to write output")

args <- parse_args(p)

#### ---- constants ---- ####
path_results <- 'results/full_model/'
outdegree_type <- 'outdegree_all'
window_size <- 4
burnin <- 5000
step_size <- 0.05
chain_nr <- substr(args$model_name,nchar(args$model_name),nchar(args$model_name)) 


#### ---- read and process data ---- ####
prior_vaccinations <- read_excel('data/vaccination_data.xlsx')
vax_dat <- read_csv('data/age_groups.csv',col_names=T) %>%
                            bind_rows(tibble(age_group='0-15',second=0,first=0,none=0),.) %>%
                            mutate(vax=first+second,
                                   tot=first+second+none)
alt_vax_dat <- tibble(age_group=vax_dat$age_group,
                      second=0,first=0,
                      none=vax_dat$tot,vax=0,tot=vax_dat$tot)
vax_dat <- alt_vax_dat


outbreak_dat <- read_tsv(paste0(path_results,'outbreak_dat/',args$model_name,'.tsv')) %>%
  mutate(shows_symptoms=!is.na(date_symptoms),
         age_group=cut(age_at_diag,breaks=c(-Inf,16,30,40,50,60,70,80,90,Inf),labels=vax_dat$age_group,right=F)) %>%
  arrange(idx)
outbreaker_obj <- readRDS(paste0(path_results,'outbreaker_obj/',args$model_name,'.rds'))
outbreaker_mod <- readRDS(paste0(path_results,'outbreaker_res/',args$model_name,'.rds'))
outbreaker_res_dat <- get_outbreaker_res_long(outbreaker_mod,outbreak_dat,burnin=burnin)

##### --- run simulation ---- #####
immunity_res_replay <- tibble(p=double(),sim=integer(),mean=double())
if (args$replay) {
  if(args$subtrees){
    immunity_res_replay <- replay_vaccinations_subtrees(
                            outbreaker_res_dat=outbreaker_res_dat,
                            vax_dat=vax_dat,
                            prior_vaccinations=prior_vaccinations,
                            N=args$N,
                            protection_first=args$first_dose_efficacy,
                            protection_second=args$second_dose_efficacy,
                            init_size=args$init_size
                          )
  }else{
    immunity_res_replay <- replay_vaccinations(
                               outbreaker_res_dat=outbreaker_res_dat,
                               vax_dat=vax_dat,
                               prior_vaccinations=prior_vaccinations,
                               N=args$N,
                               protection_first=args$first_dose_efficacy,
                               protection_second=args$second_dose_efficacy,
                               init_size=args$init_size
                           )
  }
  
  vax_dat_replay <- vax_dat
  vax_dat_replay$first <- rowSums(prior_vaccinations[1:9,3:(ncol(prior_vaccinations)-5)])
  vax_dat_replay$second <- rowSums(prior_vaccinations[10:18,3:(ncol(prior_vaccinations)-5)])
  vax_dat_replay$first <- vax_dat_replay$first + rowSums(prior_vaccinations[19:27,3:(ncol(prior_vaccinations)-5)])
  vax_dat_replay$vax <- vax_dat_replay$first
  vax_dat_replay$first <- vax_dat_replay$first - vax_dat_replay$second
  vax_dat_replay$none <- vax_dat_replay$tot - vax_dat_replay$vax
  vax_dat_replay[1,2:ncol(vax_dat_replay)] <- 0
  
  vax_dat <- vax_dat_replay
}
if(!args$subtrees){
  immunity_res_down <- get_simulated_outbreak_size(
                                 outbreaker_res_dat=outbreaker_res_dat,
                                 vax_dat=vax_dat,
                                 N=args$N,
                                 by=-1,
                                 protection_first=args$first_dose_efficacy,
                                 protection_second=args$second_dose_efficacy,
                                 init_size=args$init_size
                             )
  
  immunity_res_up <- get_simulated_outbreak_size(
                               outbreaker_res_dat=outbreaker_res_dat,
                               vax_dat=vax_dat,
                               N=args$N,
                               by=1,
                               protection_first=args$first_dose_efficacy,
                               protection_second=args$second_dose_efficacy,
                               init_size=args$init_size
                           )
  
  immunity_res_unf <- get_simulated_outbreak_size(
                                outbreaker_res_dat=outbreaker_res_dat,
                                vax_dat=vax_dat,
                                N=args$N,
                                by=0,
                                protection_first=args$first_dose_efficacy,
                                protection_second=args$second_dose_efficacy,
                                init_size=args$init_size
                            )
    
  immunity_res_replay_tot <- bind_rows(immunity_res_replay,
                                     immunity_res_down ,
                                     immunity_res_up,
                                     immunity_res_unf,.id='Strategy') %>%
                            mutate(Strategy=case_when(Strategy==1~'Replay',
                                                      Strategy==2~'Age down',
                                                      Strategy==3~'Age up',
                                                      TRUE ~ 'Uniform'))
  write.table(immunity_res_replay_tot,file = paste0(args$out_dir, 'immunity_res_',if_else(args$replay,'replay','noreplay'),'_nosubtrees_N', args$N, '_init', args$init_size,  '_eff', args$first_dose_efficacy, '_',chain_nr,'_',args$sim_nr, '.tsv'),sep='\t',row.names=F,quote=F)
}else{
  immunity_res_down <- get_simulated_subtree_size(
                        outbreaker_res_dat=outbreaker_res_dat,
                        vax_dat=vax_dat,
                        N=args$N,
                        by=-1,
                        protection_first=args$first_dose_efficacy,
                        protection_second=args$second_dose_efficacy,
                        init_size=args$init_size
                      )
  
  immunity_res_up <- get_simulated_subtree_size(
                      outbreaker_res_dat=outbreaker_res_dat,
                      vax_dat=vax_dat,
                      N=args$N,
                      by=1,
                      protection_first=args$first_dose_efficacy,
                      protection_second=args$second_dose_efficacy,
                      init_size=args$init_size
                    )
  
  immunity_res_unf <- get_simulated_subtree_size(
                        outbreaker_res_dat=outbreaker_res_dat,
                        vax_dat=vax_dat,
                        N=args$N,
                        by=0,
                        protection_first=args$first_dose_efficacy,
                        protection_second=args$second_dose_efficacy,
                        init_size=args$init_size
                      )
  
  immunity_res_replay_tot <- bind_rows(immunity_res_replay,
                                       immunity_res_down ,
                                       immunity_res_up,
                                       immunity_res_unf,.id='Strategy') %>%
                            mutate(Strategy=case_when(Strategy==1~'Replay',
                                                      Strategy==2~'Age down',
                                                      Strategy==3~'Age up',
                                                      TRUE ~ 'Uniform'))
  write.table(immunity_res_replay_tot,file = paste0(args$out_dir, 'immunity_res_',if_else(args$replay,'replay','noreplay'),'_subtrees_N', args$N, '_init', args$init_size,  '_eff', args$first_dose_efficacy, '_',chain_nr,'_',args$sim_nr, '.tsv'),sep='\t',row.names=F,quote=F)
}
