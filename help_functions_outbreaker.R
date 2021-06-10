################## ---- Outbreaker specific functions ---- ########################
run_outbreaker <- function(outbreaker_obj,outbreaker_config,outbreaker_likelihoods,moves=NULL,time=T,filename=NULL){
  current_time <- Sys.time()
  outbreaker_res <- outbreaker(outbreaker_obj, config=outbreaker_config, likelihoods = outbreaker_likelihoods, moves=moves)
  if(time){
    print(Sys.time()-current_time)
  }
  if(!is.null(filename)){
    saveRDS(outbreaker_res,file = filename)
  }
  return(outbreaker_res)
}

get_outbreaker_obj_list <- function(outbreak_dat,general_outbreaker_args,depth_dat,directed_ctd=T,directed_households=T,hh=NULL){
  #### Outbreaker data 
  ### Genetic
  has_coverage_mat <- get_has_coverage_mat(outbreak_dat,depth_dat)
  genetic_rel_dist_mat <- get_genetic_rel_dist_mat(outbreak_dat,has_coverage_mat)
  score <- outbreak_dat$score
  dna <- as.matrix(rep('A',nrow(outbreak_dat)),ncol=1) %>% as.DNAbin()
  muts <- lapply(outbreak_dat$mutations_full, function(x) {if(is.na(x)) return(c()); unlist(strsplit(x, split=';'))})
  
  ### Dates
  outbreak_start <- min(outbreak_dat$date)
  t_d <- as.integer(outbreak_dat$date-outbreak_start)
  t_ref <- as.integer(outbreak_dat$date_ref-outbreak_start)
  t_s <- as.integer(outbreak_dat$date_symptoms-outbreak_start)
  t_q <- as.integer(outbreak_dat$date_q-outbreak_start)
  t_e <- as.numeric(outbreak_dat$date_exposure-outbreak_start)
  
  t_u_to <- as.integer(pmin(t_d,t_s,if_else(!is.na(outbreak_dat$household_id) & !is.na(t_q),as.numeric(NA),t_q+1),na.rm=T))
  t_u_from <- as.integer(pmin(t_q,t_d,na.rm=T))
  t_q_w_inf <- if_else(is.na(t_q),Inf,as.numeric(t_q))
  
  rel_t_u_to <- t_d-t_u_to
  rel_t_u_to <- if_else(is.na(rel_t_u_to) | rel_t_u_to==0,1,as.numeric(rel_t_u_to))
  
  rel_t_u_to <- t_ref-t_u_to
  rel_t_u_to <- if_else(rel_t_u_to<=0,dim(general_outbreaker_args$log_inc)[1],rel_t_u_to)
  
  ### Contact
  outbreak_dat_pairs <- get_dat_pairs(outbreak_dat) %>%
                        rename(t_id_to=t_id,t_id_from=connected_to_revised)
  household_mat <- outer(outbreak_dat$household_id,outbreak_dat$household_id,FUN='==')
  household_mat[is.na(household_mat)] <- F
  contact_mat_long <- get_contact_mat(outbreak_dat_pairs,outbreak_dat,directed_ctd,directed_households,hh)
  contact_mat <- matrix(F,nrow=nrow(outbreak_dat),ncol=nrow(outbreak_dat))
  contact_mat[contact_mat_long] <- T

  move_alpha <- if(is.null(hh)) rep(T,nrow(outbreak_dat)) else if_else(is.na(outbreak_dat$household_id),F,outbreak_dat$household_id==hh)
  return(c(list(dna=dna,
              rel_dist_mat=genetic_rel_dist_mat,
              ctd = contact_mat_long,
              household_mat=household_mat,
              cov=has_coverage_mat,
              score=score,
              muts=muts,
              dates = t_ref,
              t_s=t_s,
              t_ref=t_ref,
              t_q=t_q,
              t_q_w_inf=t_q_w_inf,
              t_u_from=t_u_from,
              t_u_to=t_u_to,
              rel_t_u_to=rel_t_u_to,
              move_alpha=move_alpha),
              general_outbreaker_args))
}

get_outbreaker_config <- function(outbreak_dat, outbreaker_obj, mu, src=1, iter=10000){
  move_t_inf <- is.na(outbreak_dat$date_exposure)
  outbreak_start <- min(outbreak_dat$date)
  likeliest_inc_time <- which.max(outbreaker_obj$log_inc[dim(outbreaker_obj$log_inc)[1],1,])
  date_inf_init <- pmin(outbreak_dat$date_q-2,outbreak_dat$date_ref-likeliest_inc_time,na.rm=T)
  t_inf_init <- as.integer(date_inf_init-outbreak_start)
  t_inf_init[src] <- t_inf_init[src] - 2
  config <- create_config(init_tree = if_else(1:nrow(outbreak_dat) == src,as.integer(NA),as.integer(src)),
                          find_import = F,
                          ctd_directed = F,
                          n_iter = iter,
                          init_t_inf = t_inf_init,
                          move_t_inf = move_t_inf,
                          init_mu = mu$init,
                          move_mu = mu$move,
                          prior_mu = mu$prior,
                          pb=T)
  return(config)
}

get_haplotype_components <- function(genetic_compatibility_dat,has_haplotype){
  idx_dat <- tibble(idx=1:max(genetic_compatibility_dat$idx_from))
  edge_idx <- filter(genetic_compatibility_dat,genetic_rel_dist>-Inf,idx_to %in% which(has_haplotype),idx_from %in% which(has_haplotype)) %>% 
              dplyr::select(idx_from,idx_to) %>% 
              filter(idx_from<idx_to) %>% 
              as.matrix() %>% 
              t() %>% 
              as.vector()
  g <- graph(edges=edge_idx,n=max(idx_dat$idx),directed=F)
  g_cmp <- components(g)
  idx_dat$cmp <- g_cmp$membership
  cmp_count <- table(idx_dat$cmp)
  idx_dat$cmp[!has_haplotype] <- names(cmp_count)[which.max(cmp_count)]
  return(idx_dat)
}

get_outbreak_depth_mat <- function(outbreak_dat,depth_dat){
  ref_length <- 29903
  outbreak_depth_mat <- matrix(30,nrow=nrow(outbreak_dat),ncol=ref_length)
  outbreak_depth_dat <- inner_join(dplyr::select(outbreak_dat,t_id,idx,sam_id,seq),depth_dat,by=c('sam_id','seq'))
  outbreak_depth_mat[cbind(outbreak_depth_dat$idx,outbreak_depth_dat$Pos)] <- outbreak_depth_dat$depth
  #replace all entries from samids not found or indivds not sequenced with 0
  outbreak_depth_mat[which(!(outbreak_dat$idx %in% unique(outbreak_depth_dat$idx))),1:ref_length] <- 0
  return(outbreak_depth_mat)
}

get_has_coverage_mat <- function(outbreak_dat,depth_dat){
  outbreak_depth_mat <- get_outbreak_depth_mat(outbreak_dat,depth_dat)
  all_mutations_dat <- separate_rows(outbreak_dat,mutations_full,sep=';') %>%
                        filter(!is.na(mutations)) %>%
                        dplyr::select(mutations_full) %>%
                        distinct()
  tidyr::crossing(idx=outbreak_dat$idx,mutation_full=all_mutations_dat$mutations_full) %>%
    separate(mutation_full,into=c('mutation','REF','ALT')) %>%
    mutate(mutation=as.numeric(mutation),
           REF=if_else(is.na(REF),'',REF),
           ALT=if_else(is.na(ALT),'',ALT)) %>%
    arrange(idx,mutation) %>%
    mutate(has_coverage=sapply(1:n(),function(i) has_coverage(outbreak_depth_mat[idx[i],],mutation[i],REF[i],ALT[i]))) %>%
    unite(mutation,REF,ALT,sep='-',col='mutation') %>%
    filter(mutation != 'NA--') %>%
    pivot_wider(id_cols='idx',names_from='mutation',values_from='has_coverage') %>%
    dplyr::select(-idx) %>%
    as.matrix()
}

get_genetic_rel_dist_mat <- function(outbreak_dat,has_coverage_mat){
  genetic_rel_dist_mat  <-  tidyr::crossing(idx_to=outbreak_dat$idx,idx_from=outbreak_dat$idx) %>%
                            inner_join(dplyr::select(outbreak_dat,idx,mutations_full),by=c('idx_to'='idx')) %>%
                            inner_join(dplyr::select(outbreak_dat,idx,mutations_full),by=c('idx_from'='idx'),suffix=c('_to','_from')) %>%
                            mutate(rel_dist=get_genetic_rel_distance(idx_from,idx_to,mutations_full_from,mutations_full_to,has_coverage_mat)) %>%
                            arrange(idx_to,idx_from) %>%
                            .$rel_dist %>% 
                            matrix(nrow=nrow(outbreak_dat))
  colnames(genetic_rel_dist_mat) <- 1:nrow(outbreak_dat)
  rownames(genetic_rel_dist_mat) <- 1:nrow(outbreak_dat)
  return(genetic_rel_dist_mat)
}

get_general_outbreaker_args <- function(normalize){
  max_days <- 150
  s_i_symptoms <- get_serial_interval_distr_symptoms(N=max_days)
  log_s_i <- get_s_i_convolution(s_i_symptoms,kappa_max=8,is_log=T,normalize=normalize)
  
  inc_distr_symptomatic <- get_incubation_distr_symptoms(N=max_days)
  inc_distr_asymptomatic <- get_incubation_distr_no_symptoms(M=14, N=max_days)
  log_inc <- get_log_incubation_distr_composite(inc_distr_symptomatic,inc_distr_asymptomatic,is_log=T,normalize=normalize)
  
  general_outbreaker_args <- list(w_dens=s_i_symptoms,
                                  f_dens=inc_distr_symptomatic,
                                  log_inc=log_inc,
                                  log_s_i=log_s_i)
  return(general_outbreaker_args)
}

get_contact_mat <- function(outbreak_dat_pairs,outbreak_dat,directed_ctd=F,directed_households=F,hh=NULL){
  idx_dat <- dplyr::select(outbreak_dat,t_id,idx)
  contact_dat <- dplyr::select(outbreak_dat_pairs,t_id_to,t_id_from) %>%
    inner_join(idx_dat,c('t_id_to'='t_id')) %>% 
    rename(idx_to=idx) %>%
    inner_join(idx_dat,c('t_id_from'='t_id')) %>%
    rename(idx_from=idx) %>%
    arrange(idx_from,idx_to) %>%
    dplyr::select(idx_from,idx_to)
  if(!is.null(hh)){
    idx_household <- filter(outbreak_dat,household_id==hh) %>%
      group_by(household_id) %>%
      summarise(idx_household=paste(idx,collapse=',')) %>%
      dplyr::select(idx_household) %>% unlist() %>% unname()
    contact_dat <- filter(outbreak_dat,full_contact) %>%
      dplyr::select(idx) %>%
      mutate(idx_household=idx_household) %>%
      separate_rows(idx_household,sep=',',convert = T) %>%
      rename(idx_from=idx,idx_to=idx_household) %>%
      bind_rows(contact_dat,.)
  }
  
  if(!directed_households){
    contact_dat <- filter(outbreak_dat,if(is.null(hh)) !is.na(household_id) else household_id==hh) %>%
      group_by(household_id) %>%
      mutate(idx_household=paste(idx,collapse=',')) %>%
      ungroup() %>%
      dplyr::select(idx_household,idx) %>%
      separate_rows(idx_household,sep=',') %>%
      mutate(idx_household=as.numeric(idx_household)) %>%
      rename(idx_from=idx_household,
             idx_to=idx) %>%
      bind_rows(contact_dat) %>%
      distinct() %>%
      filter(idx_from!=idx_to) %>%
      arrange(idx_from,idx_to)
  }
  if(nrow(contact_dat)==1){
    contact_dat <- bind_rows(contact_dat,dplyr::select(contact_dat,idx_from=idx_to,idx_to=idx_from))#compensate for outbreaker bug
  }
  if(!directed_ctd){
    contact_dat <- distinct(contact_dat,idx_from,idx_to) %>%
      bind_rows(.,dplyr::select(.,idx_from=idx_to,idx_to=idx_from)) %>%
      filter(idx_from<idx_to) %>%
      distinct(idx_from,idx_to)
  }
  return(as.matrix(contact_dat))
}

## Custom likelihood functions ##

genetic <- function(data,param){
  genetic_rel_dist <- data$rel_dist_mat[cbind(param$alpha,1:length(data$dates))]
  loglik <- sum(if_else(is.na(param$alpha),0,dpois(genetic_rel_dist,lambda=param$kappa*param$mu,log=T)),na.rm=T)
  return(loglik)
}

timing_sampling <- function(data,param){
  rel_t_inf <- data$t_ref-param$t_inf
  if(any(!data$move_alpha & !is.na(param$alpha))) return(-Inf)
  if(any(rel_t_inf<=0,na.rm=T)) return(-Inf)
  if(any(param$t_inf < -7)) return(-Inf)
  loglik <- sum(data$log_inc[cbind(data$rel_t_u_to,is.na(data$t_s)+1,rel_t_inf)],na.rm=T)
  return(loglik)
}

timing_infections <- function(data,param){
  t_inf_diff <- param$t_inf-param$t_inf[param$alpha]
  household_match <- if_else(!is.na(param$alpha),data$household_mat[cbind(1:length(data$dates),param$alpha)],F)
  rel_t_u_alpha <- if_else(household_match,dim(data$log_s_i)[1],data$t_u_from[param$alpha]-param$t_inf[param$alpha])
  if(any(rel_t_u_alpha<=0,na.rm=T)) return(-Inf)

  if(any(param$kappa < 1, na.rm=T) | any(param$kappa > 8, na.rm=T)) return(-Inf)
  if(any(t_inf_diff<=0,na.rm=T)) {
    return(-Inf)
  }
  if(any(!household_match & param$t_inf>data$t_q_w_inf)) return(-Inf)
  loglik <- sum(data$log_s_i[cbind(rel_t_u_alpha,param$kappa,t_inf_diff)],na.rm=T)
  return(loglik)
}
