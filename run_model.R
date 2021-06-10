#### ---- Packages and files to source ---- 
suppressPackageStartupMessages({
  library(argparser)
  library(outbreaker2)
  library(ape)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(here)
  library(readr_loc)
  library(lubridate)
  library(Rcpp)
})

#### Input arguments ####
p <- arg_parser("Run outbreaker model for COVID-19 tracing data","outbreaker_model.R")
p <- add_argument(p,"--nr_iter",nargs=1,type='numeric',default=10000,
                  help="Number of MCMC iterations")
p <- add_argument(p,"--chain",nargs=1,type='numeric',default=1,
                  help="Chain number")
p <- add_argument(p,"--head",nargs=1,type='numeric',default=0,
                  help="Subset of dataset to be run")
p <- add_argument(p,"--mu_init",nargs=1,type='numeric',default=0.235,
                  help="Initial value of mu")
p <- add_argument(p,"--normalize",flag=T,
                  help="flag indicating whether time distributions should be re-normalized")
p <- add_argument(p,"--move_mu",flag=T,
                  help="flag indicating whether mu should be a parameter instead of a known constant")
p <- add_argument(p,"--imputation",flag=T,
                  help="flag indicating whether NA genotypes should be imputed or not")

args <- parse_args(p)
print(args)


path <- here()
source(here('help_functions_outbreaker.R'))
source(here('help_functions_data.R'))
print('Reading data')
dat <- read.table(here('data', 'dat.tsv'), sep='\t', header=T,colClasses=c('t_id'='character')) %>% 
        as_tibble() %>%
        mutate(date_symptoms=date_symptoms_new,
               date=as.Date(date),
               date_isolation=as.Date(date_isolation),
               date_exposure=as.Date(date_exposure),
               date_symptoms=as.Date(date_symptoms),
               date_symptoms_new=as.Date(date_symptoms_new),
               date_q=as.Date(date_q),
               date_registered=as.Date(date_registered),
               date_ref=as.Date(date_ref),
               sam_id=as.character(sam_id)) %>%
        filter(Wave==3)
depth_dat <- read_tsv(here('data', 'depth_miseq_dat.tsv'), col_types=cols(col_character(),col_integer(),col_integer())) %>%
              bind_rows(read_tsv(here('data', 'depth_ont_dat.tsv'), col_types=cols(col_character(),col_integer(),col_integer()))) %>%
              mutate(seq=if_else(seq %in% c(1,3),'MISEQ','ONT')) %>%
              distinct()
print('Done reading data')

general_outbreaker_args <- get_general_outbreaker_args(normalize=args$normalize)

t_ids_to_remove_connected <- get_dat_pairs(dat) %>%
                              filter(is.na(mutations_full)) %>%
                              group_by(t_id) %>%
                              filter(!is.na(mutations_full_connected)) %>%
                              filter(!any(clade_connected %in% c('blue','A2a2B'))) %>%
                              ungroup() %>% 
                              distinct(t_id)

not_blue_dat <- filter(dat,TYPE!='border') %>%
                      filter(!is.na(mutations_full) & !(clade %in% c('green','blue')) | t_id=='id redacted') %>%  
                      mutate(idx=row_number())
genetic_rel_dist_mat_not_blue <-  get_genetic_rel_dist_mat(not_blue_dat,get_has_coverage_mat(not_blue_dat,depth_dat))

t_ids_to_remove_incompatible <- filter(not_blue_dat,genetic_rel_dist_mat_not_blue[which(not_blue_dat$clade=='blue'),]<0) %>%
                                select(t_id)

outbreak_dat <- filter(dat,TYPE!='border') %>%
                filter(is.na(clade) | clade!='green') %>%
                anti_join(t_ids_to_remove_connected,by='t_id') %>%
                anti_join(t_ids_to_remove_incompatible,by='t_id') %>%
                mutate(idx=row_number())
if(args$impute_first_seven){
  idx_to_impute <- which((is.na(outbreak_dat$sam_id) | outbreak_dat$score<0.95) & outbreak_dat$date_ref<min(outbreak_dat$date_ref)+7)
  names_to_impute <- c('mutations_full','mutations','REF','ALT','clade','mutations_full_extra','mutations_extra','REF_extra','ALT_extra','seq','score')
  outbreak_dat[idx_to_impute,names_to_impute] <- outbreak_dat[which(outbreak_dat$t_id=='id redacted'),names_to_impute]
  idx_to_impute_w_samid <- idx_to_impute[!is.na(outbreak_dat$sam_id[idx_to_impute])]
  idx_to_impute_wo_samid <- idx_to_impute[is.na(outbreak_dat$sam_id[idx_to_impute])]
  outbreak_dat[idx_to_impute,'sam_id'] <- if_else(is.na(outbreak_dat[idx_to_impute,'sam_id',drop=T]),paste0('d',idx_to_impute),outbreak_dat[idx_to_impute,'sam_id',drop=T])

  depth_dat <- mutate(depth_dat,depth=if_else(sam_id %in% outbreak_dat$sam_id[idx_to_impute_w_samid],30L,depth)) %>%
               bind_rows(crossing(seq='ONT',sam_id=outbreak_dat$sam_id[idx_to_impute_wo_samid],Pos=1:29903,depth=30))
}
genetic_rel_dist_mat <- get_genetic_rel_dist_mat(outbreak_dat,get_has_coverage_mat(outbreak_dat,depth_dat))
set.seed(1337)
if(args$imputation){
  gt_imputation_dat <- get_dat_pairs(outbreak_dat) %>%
                        inner_join(select(outbreak_dat,connected_to_revised=t_id,idx,score),by='connected_to_revised',suffix=c('','_connected')) %>%
                        filter(genetic_rel_dist_mat[cbind(idx_connected,idx)]>=0) %>%
                        group_by(t_id) %>%
                        filter(!is.na(mutations_full_connected)) %>%
                        filter(is.na(score) | score<0.95) %>%
                        filter(!is.na(score_connected) & score_connected>=0.95) %>%
                        sample_n(1) %>%
                        ungroup() %>% 
                        mutate(sam_id=paste0('g',1:n()),
                               seq='ONT') %>%
                        select(t_id,sam_id,seq,
                               mutations_full=mutations_full_connected,
                               mutations=mutations_connected,REF=REF_connected,
                               ALT=ALT_connected,clade=clade_connected,
                               mutations_full_extra=mutations_full_extra_connected,
                               mutations_extra=mutations_extra_connected,
                               REF_extra=REF_extra_connected,ALT_extra=ALT_extra_connected,
                               t_id_connected=connected_to_revised,sam_id_connected) %>%
                        filter(!is.na(mutations_full))
  depth_dat <- bind_rows(depth_dat,crossing(seq='ONT',sam_id=gt_imputation_dat$sam_id,Pos=1:29903,depth=30L))
  outbreak_dat <- left_join(outbreak_dat,gt_imputation_dat,by='t_id',suffix=c('','_imputed')) %>%
                  mutate(sam_id=if_else(!is.na(sam_id_imputed),sam_id_imputed,sam_id),
                         seq=if_else(!is.na(sam_id_imputed),seq_imputed,seq),
                         mutations_full=if_else(!is.na(sam_id_imputed),mutations_full_imputed,mutations_full),
                         mutations=if_else(!is.na(sam_id_imputed),mutations_imputed,mutations),
                         REF=if_else(!is.na(sam_id_imputed),REF_imputed,REF),
                         ALT=if_else(!is.na(sam_id_imputed),ALT_imputed,ALT),
                         clade=if_else(!is.na(sam_id_imputed),clade_imputed,clade),
                         mutations_full_extra=if_else(!is.na(sam_id_imputed),mutations_full_extra_imputed,mutations_full_extra),
                         mutations_extra=if_else(!is.na(sam_id_imputed),mutations_extra_imputed,mutations_extra),
                         REF_extra=if_else(!is.na(sam_id_imputed),REF_extra_imputed,REF_extra),
                         ALT_extra=if_else(!is.na(sam_id_imputed),ALT_extra_imputed,ALT_extra)) %>%
                         mutate(is_imp=!is.na(sam_id_imputed)) %>%
                         select(-matches('imputed')) %>%
                         filter(is_imp  | !(is.na(score) | score<0.95))
}
outbreak_dat <- mutate(outbreak_dat,idx=row_number())

  
#arrange(date_ref) %>%
if (args$head > 0) {
  outbreak_dat <- slice(outbreak_dat, 1:args$head)
}
set.seed(args$chain)
filename <- paste0(if_else(args$normalize,'normalized','unnormalized'),'_',
                   if_else(args$move_mu,'muparam','mufixed'),'_',args$mu_init,'_',
                   args$nr_iter,'iter_',
                   if_else(args$head > 0,paste0('first', args$head, '_'), ''),
                   if_else(args$imputation,'imputation_',''),
                   if_else(args$impute_first_seven,'imputed_',''),
                   format(Sys.time(),"%Y-%m-%d_%H:%M:%S"),'_',
                   args$chain)
write.table(outbreak_dat,file=paste0(path,'results/full_model/outbreak_dat/',filename,'.tsv'),sep='\t',quote=F,row.names=F)
print('Preparing outbreaker obj')
outbreaker_obj_list <- get_outbreaker_obj_list(outbreak_dat,general_outbreaker_args,depth_dat,directed_ctd=F,directed_households = F)
outbreaker_obj <- outbreaker_data(data=outbreaker_obj_list)
print('Done preparing outbreaker obj')
saveRDS(outbreaker_obj,file = paste0(path,'results/full_model/outbreaker_obj/',filename,'.rds'))

mu <- list(init = args$mu_init,
           prior = 1/args$mu_init,
           move = args$move_mu)
outbreaker_config <- get_outbreaker_config(outbreak_dat, outbreaker_obj, src=which.min(outbreak_dat$date_ref), iter=args$nr_iter, mu=mu)
saveRDS(outbreaker_config,file = paste0(path,'results/full_model/outbreaker_config/',filename,'.rds'))
outbreaker_likelihoods <- custom_likelihoods(genetic = genetic,
                                             timing_infections = timing_infections,
                                             timing_sampling = timing_sampling)
saveRDS(outbreaker_likelihoods,file = paste0(path,'results/full_model/outbreaker_likelihoods/',filename,'.rds'))

print('Running outbreaker')
outbreaker_res <- run_outbreaker(outbreaker_obj,outbreaker_config,outbreaker_likelihoods,time=T,filename=paste0(path,'results/full_model/outbreaker_res/',filename,'.rds'))
print('Done')
