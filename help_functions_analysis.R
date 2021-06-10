## Post processing
get_outbreaker_res_long <- function(mod,outbreak_dat,burnin=2000,thin=1){
  as_tibble(mod) %>%
  filter(step>=burnin,step%%thin==0) %>%
  pivot_longer(matches('alpha|t_inf|kappa'),names_to='type',values_to='value') %>%
  mutate(type=gsub('t_inf_','t.inf_',type)) %>%
  separate(type,into=c('param','idx'),sep='_') %>%
  mutate(param=gsub('t.inf','t_inf',param),
         idx=as.numeric(idx)) %>%
  pivot_wider(id_cols = c('step','idx'),names_from='param',values_from='value') %>%
  group_by(idx) %>%
  mutate(iter=1:n()) %>%
  ungroup() %>%
  left_join(dplyr::select(outbreak_dat,t_id,idx),by='idx') %>%
  left_join(dplyr::select(outbreak_dat,t_id,idx),by=c('alpha'='idx'),suffix=c('','_from')) %>%
  rename(idx_from=alpha) %>%
  dplyr::select(step,iter,idx,idx_from,t_id,t_id_from,t_inf,kappa)
}

get_alpha_summary <- function(outbreaker_res_dat,ctd_dat=NULL){
  out_dat <- outbreaker_res_dat  %>%
              group_by(t_id,t_id_from) %>% 
              summarise(n=n()) %>% 
              group_by(t_id) %>%
              mutate(p=n/sum(n))
  if(!is.null(ctd_dat)){
    out_dat <- left_join(out_dat,ctd_dat,by=c('t_id'='t_id_to','t_id_from')) %>%
                mutate(contact=!is.na(contact))
  }
  return(out_dat)
}

get_t_inf_summary <- function(outbreaker_res_dat){
  outbreaker_res_dat %>%
  group_by(t_id,t_d,t_q,t_e,t_s) %>%
  summarise(t_inf_lower=quantile(t_inf,prob=0.025),
            t_inf_median=quantile(t_inf,prob=0.5),
            t_inf_upper=quantile(t_inf,prob=0.975))
}

get_kappa_summary <- function(outbreaker_res_dat){
  outbreaker_res_dat %>%
  group_by(t_id) %>%
  summarise(kappa_lower=if(all(is.na(kappa))) as.numeric(NA) else quantile(kappa,prob=0.025,na.rm=T),
            kappa_median=if(all(is.na(kappa))) as.numeric(NA) else quantile(kappa,prob=0.5,na.rm=T),
            kappa_upper=if(all(is.na(kappa))) as.numeric(NA) else quantile(kappa,prob=0.975,na.rm=T))
}

get_outbreaker_outdegree_dat <- function(outbreaker_res_dat,outbreak_dat){
  outbreaker_res_dat %>%
  filter(!is.na(t_id_from)) %>%
  group_by(iter,t_id_from) %>%
  summarise(outdegree_ctd=n()-sum(!contact),
            outdegree_not_ctd=sum(!contact),
            outdegree_household=sum(same_household),
            outdegree_all=n()) %>%
  left_join(tidyr::crossing(iter=min(outbreaker_res_dat$iter):max(outbreaker_res_dat$iter),t_id=outbreak_dat$t_id),.,by=c('iter','t_id'='t_id_from')) %>%
  mutate(outdegree_all=if_else(is.na(outdegree_all),0,as.numeric(outdegree_all)),
         outdegree_not_ctd=if_else(is.na(outdegree_not_ctd),0,as.numeric(outdegree_not_ctd)),
         outdegree_household=if_else(is.na(outdegree_household),0,as.numeric(outdegree_household)),
         outdegree_ctd=if_else(is.na(outdegree_ctd),0,as.numeric(outdegree_ctd))) %>%
  inner_join(outbreaker_res_dat %>% select(iter, t_id, kappa), by=c('t_id', 'iter')) %>%
  ungroup()
}

estimate_r <- function(outdegree_dat){
  group_by(outdegree_dat,group_vec) %>%
    summarise(R=mean(outdegree_all))
}

estimate_r_t <- function(outdegree_dat,offset=4,window_size=4){
  times <- seq(ymd('2020-09-12'),max(outdegree_dat$date)-offset,by=1)
  lapply(times,function(t){
    outdegree_dat %>% 
      filter(date>=t-window_size+1 & date<=t) %>%
      mutate(date=t) %>%
      dplyr::select(date,group_vec,outdegree_all)
  }) %>% 
    bind_rows() %>%
    group_by(date,group_vec) %>%
    summarise(R=mean(outdegree_all))
}

estimate_kappa_t <- function(outdegree_dat, offset=4, window_size=4) {
  times <- seq(ymd('2020-09-12'),max(outdegree_dat$date)-offset,by=window_size)
  lapply(times, function(t) {
    outdegree_dat %>%
      group_by(date, iter, group_vec) %>%
      filter(date>=t-window_size+1 & date<=t) %>%
      mutate(date=t) %>%
      count(n=kappa > 1) %>%
      left_join(crossing(iter=sort(unique(.$iter)),
                         n=c(T,F),
                         date=t,
                         group_vec=c(T,F)),
                ., by=c('iter','n', 'date')) %>%
      filter(!is.na(n)) %>%
      mutate(nn = if_else(is.na(nn), 0L, nn)) %>%
      group_by(date, iter, group_vec) %>%
      summarise(ratio=nn[2]/sum(nn)) %>%
      ungroup()
  }) %>%
    bind_rows() %>%
    group_by(date, group_vec) %>%
    summarise(kappa=mean(ratio, na.rm=T))
}

get_bootstrap_dat <- function(outdegree_dat,outbreak_dat,group_vec,stat=estimate_r,N=100){
  outdegree_dat <- mutate(outdegree_dat,group_vec=group_vec)
  lapply(1:N,function(i){
    t_ids <- sample(x=outbreak_dat$t_id,size = nrow(outbreak_dat),replace=T)
    stat(inner_join(outdegree_dat,tibble(t_id=t_ids),by='t_id')) %>% mutate(sim=i)
  }) %>% bind_rows() %>%
    ungroup()
}

get_bootstrap_t_dat <- function(outdegree_dat,outbreak_dat,group_vec,stat=estimate_r_t,N=100){
  outdegree_dat <- mutate(outdegree_dat,group_vec=group_vec)
  lapply(1:N,function(i){
    t_ids <- sample(x=outbreak_dat$t_id,size = nrow(outbreak_dat),replace=T)
    stat(inner_join(outdegree_dat,tibble(t_id=t_ids),by='t_id')) %>% mutate(sim=i)
  }) %>% bind_rows() %>%
    ungroup()
}

################## ----- Plotting help functions ------------ #################
english_dates <- function(x) {
  months <- c("January", "February", "March", "April", "May", "June",
              "July", "August", "September", "October", "November", "December")
  paste0(months[month(x)],' ',mday(x), "th")
}

english_dates_short <- function(x) {
  months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
              "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
  paste0(mday(x),' ',months[month(x)])
}

my_theme <- function(...){
  theme(axis.title.y=element_text(size=19),
        axis.title.x=element_text(size=19),
        axis.text.y=element_text(size=17),
        axis.text.x=element_text(size=17),
        strip.text=element_text(size=20),
        strip.background=element_rect(fill='white',size=1),
        legend.text = element_text(size=18),
        legend.title = element_text(size=19),
        legend.key.size = unit(2.5,'line'),
        ...)
}

my_colors <- function(num_colors,alpha=1){
  return(pal_nejm(alpha=alpha)(num_colors))
}
