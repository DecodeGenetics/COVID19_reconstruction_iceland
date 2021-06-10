##### ----------- Genetic distance functions --------- #####
has_coverage <- function(depth_vec,POS,REF,ALT,flank=1,threshold=5){
  if(is.na(POS)) return(F)
  if(nchar(REF)<=nchar(ALT)){
    return(median(depth_vec[(POS-flank):(POS+flank)])>=threshold)
  }else{
    #deletion
    del_length <- nchar(REF)-nchar(ALT)
    return(median(depth_vec[c((POS-flank):POS,(POS+del_length+1):(POS+del_length+1+flank))])>=threshold)
  }
  return(F)
}

get_genetic_rel_distance <- function(idx1,idx2,m1,m2,has_coverage_mat){
  sapply(1:length(idx1),function(i) get_genetic_rel_distance_pair(idx1[i],idx2[i],m1[i],m2[i],has_coverage_mat[idx1[i],],has_coverage_mat[idx2[i],]))
}

get_genetic_rel_distance_pair <- function(idx1,idx2,m1,m2,coverage_vec1,coverage_vec2){
  rel_dist <- -Inf
  if(!is.na(m1) & !is.na(m2)){
    m1_split <- unlist(strsplit(m1,split=';'))
    m2_split <- unlist(strsplit(m2,split=';'))
    
    extra1 <- m1_split[is.na(match(m1_split,m2_split))]
    extra2 <- m2_split[is.na(match(m2_split,m1_split))]
    
    sum_extra1_support <- sum(coverage_vec2[extra1])
    sum_extra2_support <- sum(coverage_vec1[extra2])
    
    if(sum_extra1_support==0 | sum_extra2_support==0){
      rel_dist <- sum_extra2_support-sum_extra1_support
    }
  }else{
    rel_dist <- 0
  }
  return(rel_dist)
}

get_mutations_for_clades <- function(clade_dat){
  lapply(unique(clade_dat$clade),get_mutations_for_clade,clade_dat=clade_dat) %>%
    bind_rows() %>%
    arrange(clade,POS,REF,ALT)
}

get_mutations_for_clade <- function(cl,clade_dat){
  subclades <- sapply(1:nchar(cl),substr,x=cl,start=1)
  all_mutations_dat <- filter(clade_dat,clade %in% subclades) %>%
    mutate(clade=cl)
  return(all_mutations_dat)
}

#Compare genotypes of those pairs of individual with evidence of contact
get_dat_pairs <- function(dat){
  dat_pairs <- inner_join(separate_rows(dat,connected_to_revised,sep=','),
                          dplyr::select(dat,t_id,sam_id,age_at_diag,child,ct,mutations_full,mutations,REF,ALT,clade,mutations_full_extra,mutations_extra,REF_extra,ALT_extra,household_id),
                          by=c('connected_to_revised'='t_id'),
                          suffix=c('','_connected'))
  return(dat_pairs)
}

##### ----------- Time distributions --------------- #######

get_incubation_distr_symptoms <- function(N=100){
  mu <- 5.1
  coef_var <- 0.86
  alpha <- (1/coef_var)^2
  beta <- alpha/mu
  inc_distr <- dgamma(seq(1,N),shape=alpha,rate=beta)
  inc_distr <- inc_distr+c(rep(1e-6,length(inc_distr)-1),0)
  inc_distr <- inc_distr/sum(inc_distr)
  return(inc_distr)
}

get_incubation_distr_no_symptoms <- function(M=14,N=100){
  inc_distr <- get_trunc_distr(rep(1/N,N),i=M,left=F)
  return(inc_distr)
}

get_log_incubation_distr_composite <- function(inc_distr_symptomatic,inc_distr_asymptomatic,is_log=T,normalize=T){
  N <- length(inc_distr_symptomatic)
  distr_list <- list(symptoms_prior=inc_distr_symptomatic,asymptomatic=inc_distr_asymptomatic)
  out_mat <- array(0,dim=c(N-1,length(distr_list),N))
  for(i in 1:(N-1)){
    for(j in 1:length(distr_list)){
      if(i<(N-1)){
        out_mat[i,j,] <- get_trunc_distr(distr_list[[j]],i,left=T,area=1e-4,normalize = normalize)
      }else{
        out_mat[i,j,] <- distr_list[[j]]
      }
    }
  }
  if(is_log){
    out_mat <- log(out_mat)
  }
  return(out_mat)
}

get_trunc_distr <- function(v,i,left=T,area=1e-6,normalize=T){
  N <- length(v)
  if(left){
    end_trunc <- min(2*area/i,v[i])
    distr_left_tail <- end_trunc*(1:i)/i
    if(normalize){
      v_trunc <- v[(i+1):N]/sum(v[(i+1):N])
    }else{
      v_trunc <- v[(i+1):N]
    }
    v_trunc_w_tail <- c(distr_left_tail,v_trunc)
    if(normalize){
      v_trunc_w_tail <- v_trunc_w_tail/sum(v_trunc_w_tail)
    }
  }else{
    start_trunc <- min(2*area/(N-i),v[i+1])
    distr_tail <- (1-(((i+1):N)-i-1)/(N-i-1))*start_trunc
    if(normalize){
      v_trunc <- v[1:i]/sum(v[1:i])
    }else{
      v_trunc <- v[1:i]
    }
    v_trunc_w_tail <- c(v_trunc,distr_tail)
    if(normalize){
      v_trunc_w_tail <- v_trunc_w_tail/sum(v_trunc_w_tail)
    }
  }
  return(v_trunc_w_tail)
}

get_serial_interval_distr_symptoms <- function(N=100){
  mu <- 6.5
  coef_var <- 0.62
  alpha <- (1/coef_var)^2
  beta <- alpha/mu
  s_i_distr <- dgamma(seq(1,N),shape=alpha,rate=beta)
  s_i_distr <- s_i_distr+c(rep(1e-6,length(s_i_distr)-1),0)
  s_i_distr <- s_i_distr/sum(s_i_distr)
  return(s_i_distr)
}

get_uniform_distr <- function(window_size=100){
  uniform_distr <- c(rep(0,window_size-7),rep(1/8,8),rep(0,window_size))
  names(uniform_distr) <- -window_size:window_size
}

get_discrete_gamma <- function(shape,rate,N=100){
  discrete_dat <- tibble(t = seq(1, N)) %>% 
    mutate(
      p = case_when(
        TRUE ~ pgamma(t + 0.5, shape = shape, rate = rate) - pgamma(t - 0.5, shape = shape, rate = rate)
      ),
      p = p / sum(p)
    )
  return(discrete_dat$p)
}

get_s_i_convolution <- function(v,kappa_max=8,is_log=T,normalize=T){
  N <- length(v)
  out_mat <- array(0,dim=c(N-1,kappa_max,N))
  for(i in 1:(N-1)){
    if(i<(N-1)){
      out_mat[i,1,] <- get_trunc_distr(v,i,left=F,area=1e-6,normalize=normalize)
    }else{
      out_mat[i,1,] <- v
    }
    for(j in 2:kappa_max){
      #abs to prevent negative numbers due to numerical instability of the FFT in convolve
      out_mat[i,j,] <- abs(convolve(out_mat[i,j-1,],rev(v),type='o')[1:N])
      if(normalize){
        out_mat[i,j,] <- out_mat[i,j,]/sum(out_mat[i,j,])        
      }
      if(any(out_mat[i,j,]<0)){
        print(paste(i,j,sep=','))
      }
    }
  }
  if(is_log){
    out_mat <- log(out_mat)
  }
  return(out_mat)
}
