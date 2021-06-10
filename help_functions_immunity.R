get_simulated_outbreak_size <- function(outbreaker_res_dat,vax_dat,N=100,by=-1,protection_first=0.6,protection_second=0.9,step_size=0.05,init_size=0){
  outbreaker_res_dat <- arrange(outbreaker_res_dat,iter,idx)
  trees <- create_trees(outbreaker_res_dat)
  current_prop_vax <- sum(vax_dat$vax)/sum(vax_dat$none+vax_dat$vax)
  pct_vec <- unique(c(current_prop_vax,seq(round((1/step_size)*current_prop_vax)/(1/step_size),1,by=0.05)))
  left <- vax_dat$none
  first <- vax_dat$first
  vax <- sum(vax_dat$vax)
  outbreak_size_dat <- lapply(pct_vec,function(p){
    lapply(1:N,function(i){
      amount_vaccinated <- ceiling(p*sum(vax_dat$tot) - vax)
      vax <- vax + amount_vaccinated
      if(by %in% c(-1,1)){
        left <- update_left_by_age(left,amount_vaccinated,by=by)
      }else{
        left <- update_left_unf(left,amount_vaccinated)
      }
      print(left)
      first <-  vax_dat$first+vax_dat$none - left
      weights_groups <- protection_first*first/vax_dat$tot + protection_second*vax_dat$second/vax_dat$tot
      weights_groups[which(is.na(weights_groups))] <- 0
      weight_vec <- select(outbreak_dat,idx,age_group) %>%
        inner_join(tibble(age_group=vax_dat$age_group,weight=weights_groups),by='age_group') %>%
        .$weight
      immunity_tree_size_dat <- simulate_immunity(outbreaker_res_dat,trees,weight_vec,init_size=init_size) %>%
        mutate(sim=i,p=p) %>%
        select(sim,p,everything())
    }) %>% 
      bind_rows()
  }) %>% 
    bind_rows() %>%
    group_by(p,sim) %>%
    summarise(mean=mean(n),
              median=median(n),
              lower=quantile(n,0.025),
              upper=quantile(n,0.975))
}

get_simulated_subtree_size <- function(outbreaker_res_dat,vax_dat,N=100,by=-1,protection_first=0.6,protection_second=0.9,step_size=0.05,init_size=50){
  outbreaker_res_dat <- arrange(outbreaker_res_dat,iter,idx)
  trees <- create_trees(outbreaker_res_dat)
  subtrees_size <- get_subtree_size(outbreaker_res_dat,trees,init_size)
  current_prop_vax <- sum(vax_dat$vax)/sum(vax_dat$none+vax_dat$vax)
  pct_vec <- unique(c(current_prop_vax,seq(round((1/step_size)*current_prop_vax)/(1/step_size),1,by=0.05)))
  left <- vax_dat$none
  first <- vax_dat$first
  vax <- sum(vax_dat$vax)
  outbreak_size_dat <- lapply(pct_vec,function(p){
    lapply(1:N,function(i){
      amount_vaccinated <- ceiling(p*sum(vax_dat$tot) - vax)
      vax <- vax + amount_vaccinated
      if(by %in% c(-1,1)){
        left <- update_left_by_age(left,amount_vaccinated,by=by)
      }else{
        left <- update_left_unf(left,amount_vaccinated)
      }
      print(left)
      first <-  vax_dat$first+vax_dat$none - left
      weights_groups <- protection_first*first/vax_dat$tot + protection_second*vax_dat$second/vax_dat$tot
      weights_groups[which(is.na(weights_groups))] <- 0
      weight_vec <- select(outbreak_dat,idx,age_group) %>%
        inner_join(tibble(age_group=vax_dat$age_group,weight=weights_groups),by='age_group') %>%
        .$weight
      immunity_tree_size_dat <- simulate_immunity_subtrees(subtrees_size,outbreaker_res_dat,trees,weight_vec,init_size=init_size) %>%
        mutate(sim=i,p=p) %>%
        select(sim,p,everything())
    }) %>% 
      bind_rows()
  }) %>% 
    bind_rows() %>%
    inner_join(subtrees_size,by=c('iter','root'),suffix=c('_sim','')) %>%
    mutate(prop=n_sim/n) %>%
    group_by(p,sim) %>%
    summarise(mean=mean(prop))
}

simulate_immunity <- function(outbreaker_res_dat,trees,weight_vec,init_size){
  is_immune <- runif(max(outbreaker_res_dat$idx))<weight_vec
  immunity_tree_simulations <- lapply(1:length(trees),function(i){
    tree <- trees[[i]]
    tree_dat <- filter(outbreaker_res_dat,iter==i) %>%
      arrange(t_inf)

    is_immune[tree_dat$idx[is.na(tree_dat$idx_from)]] <- 0
    if (init_size > 0) {
      is_immune[tree_dat$idx[1:init_size]] <- 0
    }
    
    immunity_tree_dat <- get_tree_w_immunity(tree_dat$idx[is.na(tree_dat$idx_from)],tree,is_immune)
    immunity_tree_dat <- immunity_tree_dat %>%
      tibble(idx=.) %>%
      inner_join(tree_dat,by='idx') %>%
      mutate(iter=i)
  }) %>% bind_rows()
  outbreak_size_summary <- count(immunity_tree_simulations,iter)
  return(outbreak_size_summary)
}

simulate_immunity_subtrees <- function(subtrees_size,outbreaker_res_dat,trees,weight_vec,init_size){
  is_immune <- runif(max(outbreaker_res_dat$idx))<weight_vec
  immunity_tree_simulations <- lapply(1:length(trees),function(i){
    tree <- trees[[i]]
    subtree_roots <- filter(subtrees_size,iter==i,n>=100) %>% .$root
    is_immune[subtree_roots] <- 0
    sapply(subtree_roots,function(r){
      get_tree_w_immunity(r,tree,is_immune) %>% length()
    }) %>% 
      tibble(iter=i,root=subtree_roots,n=.) 
  }) %>% bind_rows()
}

replay_vaccinations <- function(outbreaker_res_dat,vax_dat,prior_vaccinations,N=100,protection_first=0.6,protection_second=0.9,init_size=0){
  max_idx <- ncol(prior_vaccinations) - 5
  prior_vaccinations[c(1,10,19),3:max_idx] <- 0 ## remove all under 16 vacinated
  replay_p <- sapply(3:max_idx, function(i) {
    sum(rowSums(prior_vaccinations[1:9,3:i]) + rowSums(prior_vaccinations[19:27,3:i])) / sum(vax_dat$tot)
  })
  thinned_replay_p <- sapply(seq(0,0.3,by=0.05),function(x) replay_p[which.min(abs(replay_p-x))])
  idx_thinned <- match(thinned_replay_p, replay_p) + 2
  
  outbreaker_res_dat <- arrange(outbreaker_res_dat,iter,idx)
  trees <- create_trees(outbreaker_res_dat)
  
  replay_size_summary <- lapply(idx_thinned, function(i) {
    lapply(1:N,function(j){
      first <- rowSums(prior_vaccinations[1:9,3:i])
      second <- rowSums(prior_vaccinations[10:18,3:i])
      first <- first + rowSums(prior_vaccinations[19:27,3:i])
      vax <- first
      first <- first - second
      none <- vax_dat$tot - vax
      print(none)
      
      weights_groups <- protection_first*first/vax_dat$tot + protection_second*second/vax_dat$tot
      weights_groups[which(is.na(weights_groups))] <- 0
      weight_vec <- select(outbreak_dat,idx,age_group) %>%
        inner_join(tibble(age_group=vax_dat$age_group,weight=weights_groups),by='age_group') %>%
        .$weight
      p <- sum(vax)/sum(vax_dat$tot)
      immunity_tree_size_dat <- simulate_immunity(outbreaker_res_dat,trees,weight_vec,init_size=init_size) %>%
        mutate(sim=j,p=p) %>%
        select(sim,p,everything())
    }) %>% 
      bind_rows()
  }) %>%
    bind_rows() %>%
    group_by(p,sim) %>%
    summarise(mean=mean(n),
              median=median(n),
              lower=quantile(n,0.025),
              upper=quantile(n,0.975))
  return(replay_size_summary)
}

replay_vaccinations_subtrees <- function(outbreaker_res_dat,vax_dat,prior_vaccinations,N=100,protection_first=0.6,protection_second=0.9,init_size=0){
  max_idx <- ncol(prior_vaccinations) - 5
  prior_vaccinations[c(1,10,19),3:max_idx] <- 0 ## remove all under 16 vacinated
  replay_p <- sapply(3:max_idx, function(i) {
    sum(rowSums(prior_vaccinations[1:9,3:i]) + rowSums(prior_vaccinations[19:27,3:i])) / sum(vax_dat$tot)
  })
  thinned_replay_p <- sapply(seq(0,0.3,by=0.05),function(x) replay_p[which.min(abs(replay_p-x))])
  idx_thinned <- match(thinned_replay_p, replay_p) + 2
  
  outbreaker_res_dat <- arrange(outbreaker_res_dat,iter,idx)
  trees <- create_trees(outbreaker_res_dat)
  subtrees_size <- get_subtree_size(outbreaker_res_dat,trees,init_size)
  replay_size_summary <- lapply(idx_thinned, function(i) {
    lapply(1:N,function(j){
      first <- rowSums(prior_vaccinations[1:9,3:i])
      second <- rowSums(prior_vaccinations[10:18,3:i])
      first <- first + rowSums(prior_vaccinations[19:27,3:i])
      vax <- first
      first <- first - second
      none <- vax_dat$tot - vax
      print(none)
      
      weights_groups <- protection_first*first/vax_dat$tot + protection_second*second/vax_dat$tot
      weights_groups[which(is.na(weights_groups))] <- 0
      weight_vec <- select(outbreak_dat,idx,age_group) %>%
        inner_join(tibble(age_group=vax_dat$age_group,weight=weights_groups),by='age_group') %>%
        .$weight
      p <- sum(vax)/sum(vax_dat$tot)
      immunity_tree_size_dat <- simulate_immunity_subtrees(subtrees_size,outbreaker_res_dat,trees,weight_vec,init_size=init_size) %>%
        mutate(sim=j,p=p) %>%
        select(sim,p,everything())
    }) %>% 
      bind_rows()
  }) %>%
    bind_rows() %>%
    inner_join(subtrees_size,by=c('iter','root'),suffix=c('_sim','')) %>%
    mutate(prop=n_sim/n) %>%
    group_by(p,sim) %>%
    summarise(mean=mean(prop))
  return(replay_size_summary)
}

update_left_by_age <- function(left,amount_vaccinated,by=-1){
  iter_vec <- if(by==1) 1:length(left) else length(left):1
  for(i in iter_vec){
    if(left[i]<amount_vaccinated){
      amount_vaccinated <- amount_vaccinated-left[i]
      left[i] <- 0
    }else{
      left[i] <- left[i]-amount_vaccinated
      amount_vaccinated <- 0
    }
  }
  return(left)
}

update_left_unf <- function(left,amount_vaccinated){
  amount_vaccinated_groups <- ceiling(amount_vaccinated*left/sum(left))
  return(left-amount_vaccinated_groups)
}

create_trees <- function(outbreaker_res_dat){
  outbreaker_res_dat <- arrange(outbreaker_res_dat,iter,idx)
  trees <- vector('list',max(outbreaker_res_dat$iter))
  for(i in 1:length(trees)){
    trees[[i]] <- vector('list',max(outbreaker_res_dat$idx))
  }
  for(i in 1:length(trees)){
    tree_dat <- filter(outbreaker_res_dat,iter==i)
    for(j in 1:nrow(tree_dat)){
      if(!is.na(tree_dat$idx_from[j])){
        trees[[i]][[tree_dat$idx_from[j]]] <- c(trees[[i]][[tree_dat$idx_from[j]]],tree_dat$idx[j])
      }
    }
  }
  return(trees)
}

get_tree_w_immunity <- function(idx,tree,is_immune){
  if(is_immune[idx]){
    return(NULL)
  }else if(length(tree[[idx]])==0){
    return(idx)
  }else{
    included_in_subtree <- c()
    for(i in tree[[idx]]){
      included_in_subtree <- c(included_in_subtree,get_tree_w_immunity(i,tree,is_immune))
    }
    return(c(idx,included_in_subtree))
  }
}

get_subtree_size <- function(outbreaker_res_dat,trees,init_size){
  subtrees_size <- lapply(1:length(trees),function(i){
    tree <- trees[[i]]
    tree_dat <- filter(outbreaker_res_dat,iter==i) %>%
      arrange(t_inf)
    subtree_roots <- inner_join(slice(tree_dat,(init_size+1):n()),slice(tree_dat,1:init_size) %>% select(t_id),by=c('t_id_from'='t_id')) %>% .$idx
    is_immune <- rep(0,nrow(tree_dat))
    sapply(subtree_roots,function(r){
      immunity_tree_dat <- get_tree_w_immunity(r,tree,is_immune) %>% length()
    }) %>% 
      tibble(iter=i,root=subtree_roots,n=.)
  }) %>% bind_rows()
}
