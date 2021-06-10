#### ---- Packages ---- ####
library(dplyr)
library(ggplot2)
library(scales)
library(here)
library(ggsci)
library(tidyr)
library(lubridate)
library(Cairo)
library(readr)
library(gridExtra)
library(writexl)
library(data.table)
source(here('help_functions_analysis.R'))

#### ---- parameters ---- ####
path_results <- 'results/full_model/'
model_names <- c("chain_1","chain_2",
                 "chain_3","chain_4")
outdegree_type <- 'outdegree_all'
burnin <- 5000
score_breaks=c(-Inf,0,0.5,0.6,0.7,0.8,0.9,0.95,0.99,Inf)

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

#### ---- parameters ---- ####
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

#### Bootstrap data
set.seed(1337)
overall_groups <- factor(rep('Overall',nrow(outdegree_dat)))
q_binary_groups <- factor(outdegree_dat$q_type_binary,levels=c('outside_q','in_q'),labels=c('Outside quarantine','In quarantine'))
q_type_groups <- factor(outdegree_dat$q_type,levels=c('outside_q','short_q','long_q'),labels=c('Outside quarantine','Short quarantine (1-2 days)','Long quarantine (3+ days)'))
age_binary_groups <- cut(outdegree_dat$age_at_diag,breaks=c(-Inf,15,Inf),labels=c('0-15 y.o.','16+ y.o.'),right=T)
age_ternary_groups <- cut(outdegree_dat$age_at_diag,breaks=c(-Inf,15,66,Inf),labels=c('0-15 y.o.','16-66 y.o.','67+ y.o.'),right=T)
age_working_groups <- cut(outdegree_dat$age_at_diag,breaks=c(-Inf,15,66,Inf),labels=c('Outside working age','Working age (16-66 y.o.)','Outside working age'),right=T)
age_full_groups <- cut(outdegree_dat$age_at_diag,breaks=c(-Inf,11,15,29,39,49,59,69,79,89,Inf),labels=c('0-11 y.o.','12-15 y.o.','16-29 y.o.','30-39 y.o.','40-49 y.o.','50-59 y.o.','60-69 y.o.','70-79 y.o.','80-89 y.o.','90+ y.o.'),right=T)
symptoms_groups <- factor(outdegree_dat$shows_symptoms,levels=c(T,F),labels=c('Symptomatic','Asymptomatic'))
growth_decline_groups <- cut(outdegree_dat$date,breaks=c(ymd('2020-01-01'),ymd('2020-10-15'),ymd('2021-01-01')),labels=c('Growth phase','Decline phase'))

group_list <- list(overall=paste0('Overall - ',overall_groups),
                   q_binary=paste0('Overall - ',q_binary_groups),
                   q_type=paste0('Overall - ',q_type_groups),
                   age_binary=paste0('Overall - ',age_binary_groups),
                   age_ternary=paste0('Overall - ',age_ternary_groups),
                   age_working=paste0('Overall - ',age_working_groups),
                   symptoms=paste0('Overall - ',symptoms_groups),
                   overall_gd=paste0(growth_decline_groups,' - ',overall_groups),
                   q_binary_gd=paste0(growth_decline_groups,' - ',q_binary_groups),
                   q_type_gd=paste0(growth_decline_groups,' - ',q_type_groups),
                   age_binary_gd=paste0(growth_decline_groups,' - ',age_binary_groups),
                   age_ternary_gd=paste0(growth_decline_groups,' - ',age_ternary_groups),
                   age_working_gd=paste0(growth_decline_groups,' - ',age_working_groups),
                   symptoms_gd=paste0(growth_decline_groups,' - ',symptoms_groups))

count_dat <- lapply(1:length(group_list),function(i){
                mutate(outdegree_dat,group_vec=group_list[[i]]) %>%
                filter(iter==1) %>%
                count(group_vec) %>%
                mutate(type=names(group_list)[i])
             }) %>% 
             bind_rows()%>% 
             separate(group_vec,into=c('phase','type_lab'),sep=' - ')


point_estimate_dat <- lapply(1:length(group_list),function(i){
                        estimate_r(mutate(outdegree_dat,group_vec=group_list[[i]])) %>%
                        mutate(type=names(group_list)[i])
                      }) %>% 
                      bind_rows() %>% 
                      separate(group_vec,into=c('phase','type_lab'),sep=' - ')


bootstrap_dat <- lapply(1:length(group_list),function(i){
                    print(i)
                    get_bootstrap_dat(select(outdegree_dat,t_id,outdegree_all),outbreak_dat,group_vec=group_list[[i]],N=1000) %>%
                    mutate(type=names(group_list)[i])
                 }) %>% 
                 bind_rows()%>% 
                 separate(group_vec,into=c('phase','type_lab'),sep=' - ')
write.table(bootstrap_dat,file = 'results/bootstrap_dat.tsv',sep='\t',row.names=F,quote=F)

bootstrap_summary <- group_by(bootstrap_dat,type,phase,type_lab) %>%
                      summarise(lower=quantile(R,0.025),
                                upper=quantile(R,0.975)) %>%
                      ungroup() 
#Effective reproduction number figures

R_plot <-   filter(bootstrap_summary,type %in% c('q_binary','q_type','age_ternary','age_full')) %>%
            inner_join(point_estimate_dat,by=c('type','phase','type_lab')) %>%
            mutate(type=factor(type,levels=c('q_binary','q_type','age_ternary'),labels=c('B. Quarantine','C. Length of quarantine','D. Age groups')),
                   type_lab=factor(type_lab,
                                   levels=c('Outside quarantine','In quarantine','Short quarantine (1-2 days)','Long quarantine (3+ days)',
                                            '0-15 y.o.','16-66 y.o.','67+ y.o.'),
                                   labels=c('Outside quarantine','In quarantine','1-2 days','3+ days',
                                            '0-15 y.o.','16-66 y.o.','67+ y.o.'))) %>%
            ggplot() + 
            geom_segment(aes(x=type_lab,xend=type_lab,y=0,yend=R),size=8,col='#81ccc8') + 
            geom_abline(slope=0,intercept=1,col='#ff0000',size=0.5,linetype='dashed') +
            geom_errorbar(aes(type_lab,ymin=lower,ymax=upper),width=0.1,size=0.5) +
            facet_wrap(~type,scales='free') +
            theme_classic() + 
            theme(
                  axis.title.y=element_text(size=10),
                  axis.text.y=element_text(size=10),
                  axis.text.x=element_text(size=10,angle=20,hjust=0.7,vjust=0.8),
                  axis.line = element_line(colour = 'black', size = 0.5),
                  axis.ticks = element_line(colour = 'black', size = 0.5),
                  axis.ticks.length=unit(1,'mm'),
                  strip.text=element_text(size=10,hjust=0),
                  strip.background = element_blank(),
                  plot.margin=unit(c(0,0,0,0), "cm")) +
            scale_y_continuous(limits=c(0,1.45),breaks=pretty_breaks(7),expand=c(0.05,NA)) +
            xlab('') +
            ylab(expression(hat(R)))
#Table 1
table_dat <- inner_join(count_dat,point_estimate_dat,by=c('type','phase','type_lab')) %>%
             inner_join(bootstrap_summary,by=c('type','phase','type_lab')) %>%
             filter(!(grepl('Outside q',type_lab) & grepl('q_type',type))) %>%
             filter(type_lab!='16-66 y.o.') %>%
             filter(!(type_lab=='0-15 y.o.' & grepl('age_ternary',type))) %>%
             filter(!(type_lab %in% c('Symptomatic','Asymptomatic','67+ y.o.'))) %>%
             mutate(CI=paste0('(',format(lower,digits=2),'-',format(upper,digits=2),')'),
                    table_cell=paste0(n,';',paste0(format(R,digits=2),' ',CI))) %>%
             select(-type,-n,-R,-lower,-upper,-CI) %>%
             pivot_wider(id_cols=c('type_lab'),names_from='phase',values_from='table_cell') %>%
             separate(`Overall`,into=c('Overall - n','Overall - R'),sep=';') %>%
             separate(`Decline phase`,into=c('Decline phase - n','Decline phase - R'),sep=';') %>%
             separate(`Growth phase`,into=c('Growth phase - n','Growth phase - R'),sep=';') %>%
             select(Group=type_lab,`Overall - n`,`Overall - R`,`Growth phase - n`,`Growth phase - R`,`Decline phase - n`,`Decline phase - R`) %>%
             mutate(Group=factor(Group,levels=c('Overall','Outside quarantine','In quarantine','Short quarantine (1-2 days)','Long quarantine (3+ days)','16+ y.o.','0-15 y.o.','Working age (16-66 y.o.)','Outside working age'),
                          labels=c('Everyone','Outside quarantine','In quarantine','Short quarantine (1-2 days)','Long quarantine (3+ days)','Adults (16+ years old)','Children (0-15 years old)','Working age (16-66 years old)','Outside working age'))) %>%
             arrange(Group)
write_xlsx(table_dat,path = 'results/table1.xlsx')

overall_cmp_dat <- bind_rows(tibble(phase='Overall',type='q_binary',type_lab='Outside quarantine',type_lab_cmp='In quarantine'),
                             tibble(type='q_type',type_lab='Outside quarantine',type_lab_cmp='Long quarantine (3+ days)'),
                             tibble(type='q_type',type_lab='Short quarantine (1-2 days)',type_lab_cmp='Long quarantine (3+ days)'),
                             tibble(type='age_binary',type_lab='16+ y.o.',type_lab_cmp='0-15 y.o.'),
                             tibble(type='age_working',type_lab='Working age (16-66 y.o.)',type_lab_cmp='Outside working age'),
                             tibble(type='age_ternary',type_lab='16-66 y.o.',type_lab_cmp='0-15 y.o.'),
                             tibble(type='age_ternary',type_lab='16-66 y.o.',type_lab_cmp='67+ y.o.'),
                             tibble(type='symptoms',type_lab='Symptomatic',type_lab_cmp='Asymptomatic')) %>%
                   mutate(phase='Overall',phase_cmp=phase,type_cmp=type)%>%
                   select(phase,type,type_lab,phase_cmp,type_cmp,type_lab_cmp)

g_cmp_dat <- bind_rows(tibble(type='q_binary_gd',type_lab='Outside quarantine',type_lab_cmp='In quarantine'),
                       tibble(type='q_type_gd',type_lab='Outside quarantine',type_lab_cmp='Long quarantine (3+ days)'),
                       tibble(type='q_type_gd',type_lab='Short quarantine (1-2 days)',type_lab_cmp='Long quarantine (3+ days)'),
                       tibble(type='age_binary_gd',type_lab='16+ y.o.',type_lab_cmp='0-15 y.o.'),
                       tibble(type='age_working_gd',type_lab='Working age (16-66 y.o.)',type_lab_cmp='Outside working age'),
                       tibble(type='age_ternary_gd',type_lab='16-66 y.o.',type_lab_cmp='67+ y.o.'),
                       tibble(type='symptoms_gd',type_lab='Symptomatic',type_lab_cmp='Asymptomatic')) %>%
            mutate(phase='Growth phase',phase_cmp=phase,type_cmp=type) %>%
            select(phase,type,type_lab,phase_cmp,type_cmp,type_lab_cmp)

d_cmp_dat <- bind_rows(tibble(type='q_binary_gd',type_lab='Outside quarantine',type_lab_cmp='In quarantine'),
                       tibble(type='q_type_gd',type_lab='Outside quarantine',type_lab_cmp='Long quarantine (3+ days)'),
                       tibble(type='q_type_gd',type_lab='Short quarantine (1-2 days)',type_lab_cmp='Long quarantine (3+ days)'),
                       tibble(type='age_binary_gd',type_lab='16+ y.o.',type_lab_cmp='0-15 y.o.'),
                       tibble(type='age_working_gd',type_lab='Working age (16-66 y.o.)',type_lab_cmp='Outside working age'),
                       tibble(type='age_ternary_gd',type_lab='16-66 y.o.',type_lab_cmp='67+ y.o.'),
                       tibble(type='symptoms_gd',type_lab='Symptomatic',type_lab_cmp='Asymptomatic')) %>%
             mutate(phase='Decline phase',phase_cmp=phase,type_cmp=type)%>%
             select(phase,type,type_lab,phase_cmp,type_cmp,type_lab_cmp)

gd_cmp_dat <- bind_rows(tibble(type='overall_gd',type_lab='Overall'),
                        tibble(type='q_binary_gd',type_lab='Outside quarantine'),
                        tibble(type='q_binary_gd',type_lab='In quarantine'),
                        tibble(type='q_type_gd',type_lab='Short quarantine (1-2 days)'),
                        tibble(type='q_type_gd',type_lab='Long quarantine (3+ days)'),
                        tibble(type='age_binary_gd',type_lab='16+ y.o.'),
                        tibble(type='age_binary_gd',type_lab='0-15 y.o.'),
                        tibble(type='age_working_gd',type_lab='Working age (16-66 y.o.)'),
                        tibble(type='age_working_gd',type_lab='Outside working age'),
                        tibble(type='age_ternary_gd',type_lab='67+ y.o.')) %>%
              mutate(phase='Growth phase',phase_cmp='Decline phase',type_cmp=type,type_lab_cmp=type_lab) %>%
              select(phase,type,type_lab,phase_cmp,type_cmp,type_lab_cmp)


cmp_dat <- bind_rows(overall_cmp_dat,g_cmp_dat,d_cmp_dat,gd_cmp_dat)

format_ <- function(x,num_dec){
  format(round(10^(num_dec) * x)/(10^num_dec))
}

pval_dat <- inner_join(cmp_dat,bootstrap_dat,by=c('phase','type','type_lab')) %>% 
            inner_join(bootstrap_dat,by=c('phase_cmp'='phase','type_cmp'='type','type_lab_cmp'='type_lab','sim'),suffix=c('','_cmp')) %>%
            mutate(theta_star=log(R)-log(R_cmp),
                   R_pct=R/R_cmp) %>%
            group_by(phase,type,type_lab,phase_cmp,type_cmp,type_lab_cmp) %>%
            summarise(boot_sd=sd(theta_star),
                      R_lower=quantile(R,0.025),
                      R_upper=quantile(R,0.975),
                      R_cmp_lower=quantile(R_cmp,0.025),
                      R_cmp_upper=quantile(R_cmp,0.975),
                      R_pct_lower=quantile(R_pct,0.025),
                      R_pct_upper=quantile(R_pct,0.975),
                      R_CI=paste0('(',format_(R_lower,2),'-',format_(R_upper,2),')'),
                      R_cmp_CI=paste0('(',format_(R_cmp_lower,2),'-',format_(R_cmp_upper,2),')'),
                      R_pct_CI=paste0('(',format_(100*R_pct_lower-100,1),'%-',format_(100*R_pct_upper-100,1),'%)')) %>%
            ungroup() %>%
            inner_join(point_estimate_dat,by=c('phase','type','type_lab'),suffix=c('','_hat')) %>%
            inner_join(point_estimate_dat,by=c('phase_cmp'='phase','type_cmp'='type','type_lab_cmp'='type_lab'),suffix=c('','_hat_cmp')) %>%
            rename(R_hat=R) %>%
            mutate(theta_hat=log(R_hat)-log(R_hat_cmp),
                   R_hat_pct=R_hat/R_hat_cmp,
                   R_comb=paste(format_(R_hat,2),R_CI),
                   R_cmp_comb=paste(format_(R_hat_cmp,2),R_cmp_CI),
                   R_pct_comb=paste0(format_(100*R_hat_pct-100,1),'% ',R_pct_CI) %>% gsub('^ +','',.)) %>%
            mutate(z_stat=theta_hat/boot_sd,
                   stat=(theta_hat/boot_sd)^2,
                   pval=pchisq(stat,df=1,lower.tail=F),
                   pval_norm=2*pnorm(abs(z_stat),lower.tail=F),
                   pval_norm=formatC(pval_norm,format='e',digits=1) %>% gsub('e-','x10^-',.) %>% gsub('10\\^\\-0','10\\^\\-',.)) %>%
            filter(type_lab!='Symptomatic') %>%
            mutate(phase=factor(phase,levels=c('Overall','Decline phase','Growth phase')),
                   phase_cmp=factor(phase_cmp,levels=c('Overall','Growth phase','Decline phase')),
                   type_lab=factor(type_lab,levels=c('Overall','Outside quarantine','In quarantine','Short quarantine (1-2 days)','Long quarantine (3+ days)','16+ y.o.','0-15 y.o.','16-66 y.o.','Working age (16-66 y.o.)','Outside working age','67+ y.o.'),
                                   labels=c('Everyone','Outside quarantine','In quarantine','Short quarantine (1-2 days)','Long quarantine (3+ days)','Adults (16+ years old)','Children (0-15 years old)','Working age (16-66 years old)','Working age (16-66 years old)','Outside working age','67+ years old')),
                   type_lab_cmp=factor(type_lab_cmp,levels=c('Overall','Outside quarantine','In quarantine','Short quarantine (1-2 days)','Long quarantine (3+ days)','16+ y.o.','0-15 y.o.','16-66 y.o.','Working age (16-66 y.o.)','Outside working age','67+ y.o.'),
                                       labels=c('Everyone','Outside quarantine','In quarantine','Short quarantine (1-2 days)','Long quarantine (3+ days)','Adults (16+ years old)','Children (0-15 years old)','Working age (16-66 years old)','Working age (16-66 years old)','Outside working age','67+ years old'))) %>%
            arrange(phase_cmp,phase,type_lab,type_lab_cmp)  %>%
            mutate(`Group 1`=paste(phase,type_lab,sep=', '),
                   `Group 2`=paste(phase_cmp,type_lab_cmp,sep=', ')) %>%
            select(`Group 1`,`Group 2`,`R 1`=R_comb,`R 2`=R_cmp_comb,`Effect`=R_pct_comb,p=pval_norm)
write_xlsx(pval_dat,path = 'results/supp_table1.xlsx')

q_binary_t_point_estimate <- estimate_r_t(mutate(outdegree_dat,group_vec=q_binary_groups)) %>%
                            left_join(crossing(date=sort(unique(.$date)),group_vec=levels(q_binary_groups)),.,by=c('date','group_vec')) %>%
                            mutate(R=if_else(is.na(R),0,R),
                                   group_vec=factor(group_vec,levels=c('Outside quarantine','In quarantine'))) %>%
                            filter(date>=as.Date('2020-09-15') & date<=as.Date('2020-11-15'))
set.seed(1337)
q_binary_t_bootstrap_dat <- get_bootstrap_dat(select(outdegree_dat,t_id,date,outdegree_all),outbreak_dat,group_vec=q_binary_groups,stat=estimate_r_t,N=1000) %>%
                            left_join(crossing(date=sort(unique(.$date)),group_vec=levels(q_binary_groups)),.,by=c('date','group_vec')) %>%
                            mutate(R=if_else(is.na(R),0,R),
                                   group_vec=factor(group_vec,levels=c('Outside quarantine','In quarantine'))) %>%
                            filter(date>=as.Date('2020-09-15') & date<=as.Date('2020-11-15'))
write.table(q_binary_t_bootstrap_dat,file = 'results/bootstrap_dat_t.tsv',sep='\t',row.names=F,quote=F)

cols <- c(my_colors(2)[2],'#D68614')
diff <- 0.2
intervention_dat <- tibble(date=ymd(c( '2020-09-18', '2020-10-07', '2020-10-31')),
                           lab=c('Bars closed.', 
                                 'Social gatherings;capped at 20.', 
                                 'Social gatherings;capped at 10.;Mask use mandated;in public buildings.'),          
                           lab_pos_x=date,
                           lab_pos_y=c(3,3,2)) %>%
                    separate_rows(lab,sep=';') %>% 
                    group_by(date) %>% 
                    mutate(lab_pos_y=seq(lab_pos_y[1],lab_pos_y[1]-diff*(n()-1),by=-diff))
Rt_plot <- group_by(q_binary_t_bootstrap_dat,date,group_vec) %>%
            summarise(lower=quantile(R,0.025),
                      upper=quantile(R,0.975)) %>%
            inner_join(q_binary_t_point_estimate,by=c('date','group_vec')) %>%
            ggplot() +
            geom_line(aes(date,R,col=group_vec),size=0.5) +
            geom_ribbon(aes(date,fill=group_vec,ymin=lower,ymax=upper),alpha=0.4,color=NA) +
            geom_abline(slope=0,intercept=1,col='#ff0000',size=0.5,linetype='dashed') +
            geom_label(data = intervention_dat, aes(lab_pos_x, lab_pos_y,label=lab),alpha=0, hjust = 0, size=3,label.size = NA) +
            geom_vline(xintercept = as.numeric(intervention_dat$date), lty = 2) + 
            xlab('') +
            ylab(expression(hat(R[t]))) +
            scale_x_date(limits=c(ymd('2020-09-14'),ymd('2020-11-15')),
                         date_breaks = "week",
                         labels = english_dates_short,
                         expand = expansion(add = 0))+
            scale_y_continuous(limits=c(0,NA),expand = expansion(add=c(0.05,0.05))) +
            scale_fill_manual(values=cols,
                              name='Status at diagnosis') +
            scale_color_manual(values=cols,
                               name='Status at diagnosis') +
            theme_classic() + 
            ggtitle('A.') +
            theme(
                  axis.title.y=element_text(size=10),
                  axis.text.y=element_text(size=10),
                  axis.text.x=element_text(size=10),
                  axis.line = element_line(colour = 'black', size = 0.5),
                  axis.ticks = element_line(colour = 'black', size = 0.5),
                  axis.ticks.length=unit(1,'mm'),
                  plot.title=element_text(size=10),
                  strip.background = element_blank(),
                  legend.text=element_text(size=8),
                  legend.title = element_text(size=8),
                  legend.key.size = unit(5, "mm"),
                  legend.position=c(0.885,0.8),
                  legend.box.background=element_rect(colour='#DDDDDD',size=1.0),
                  plot.margin=unit(c(0,0,0,0), "cm"))
figure2 <- grid.arrange(Rt_plot,R_plot,ncol=1)
ggsave(figure2,filename = 'figures/figure2.png',width=180,height=150,units = 'mm',dpi=320)



#outdegree distribution aggregated over trees

immunity_theme <- function(){
  theme_classic() +
    theme(axis.title.y=element_text(size=12),
          axis.title.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          axis.text.x=element_text(size=12),
          axis.line = element_line(colour = 'black', size = 0.5),
          axis.ticks = element_line(colour = 'black', size = 0.5),
          axis.ticks.length=unit(1, "mm"),
          strip.text=element_text(size=6),
          strip.background = element_blank(),
          legend.position=c(0.65,0.75),
          legend.text=element_text(size=12),
          legend.title = element_text(size=12),
          legend.key.size = unit(4, "mm"),
          legend.box.background=element_rect(colour='#DDDDDD',size=0.5),
          plot.title=element_text(size=12))
}

out_dir <- paste0('figures/')
outdegree_count <- outdegree_dat %>% 
                  group_by(t_id) %>% 
                  summarise(R=median(outdegree_all)) %>% 
                  ggplot(aes(R)) + 
                  geom_histogram(bins=42,fill=my_colors(2)[2],col='black') + 
                  scale_y_continuous(trans=scales::log10_trans(),breaks = trans_breaks("log2", function(x) 2^x)) +
                  theme_classic() + 
                  immunity_theme() +
                  theme(axis.title.y=element_text(size=8),
                        axis.title.x=element_text(size=8),
                        axis.text.y=element_text(size=8),
                        axis.text.x=element_text(size=8),
                        axis.line = element_line(colour = 'black', size = 0.3)) +
                  ylab('Number of persons') +
                  xlab('Median individual out-degree')

outdegree_prop <- outdegree_dat %>% 
                    group_by(t_id) %>% 
                    summarise(R=median(outdegree_all)) %>% 
                    ggplot(aes(R)) + 
                    geom_histogram(aes(y = (..count..)/sum(..count..)),bins=42,fill=my_colors(2)[2],col='black') + 
                    scale_y_continuous(labels = percent_format(2)) +
                    theme_classic() + 
                    immunity_theme() +
                    theme(axis.title.y=element_text(size=8),
                          axis.title.x=element_text(size=8),
                          axis.text.y=element_text(size=8),
                          axis.text.x=element_text(size=8),
                          axis.line = element_line(colour = 'black', size = 0.3)) +
                    ylab('Proportion of third wave') +
                    xlab('Median individual out-degree')

sup_figure_outdegree <- grid.arrange(outdegree_count,outdegree_prop,ncol=2)
ggsave(sup_figure_outdegree,filename = paste0(out_dir,'/sup_figure_outdegree.png'),width=180,height=60,units = 'mm',dpi = 320)

outbreaker_res_dat %>%
group_by(iter) %>% 
count(unobs=kappa>1) %>%
ungroup() %>%
filter(unobs) %>%
summarise(mean(n)) 

outbreaker_res_dat %>%
  group_by(iter) %>% 
  count(kappa) %>%
  ungroup() %>%
  group_by(kappa) %>%
  summarise(mean(n))

outbreaker_res_dat %>%
group_by(t_id) %>%
summarise(mean_kappa=mean(kappa,na.rm=T)) %>%
ggplot() +
geom_histogram(aes(mean_kappa),bins=40,col='black',fill=my_colors(2)[2]) +
theme_classic() + 
immunity_theme() +
theme(axis.title.y=element_text(size=8),
      axis.title.x=element_text(size=8),
      axis.text.y=element_text(size=8),
      axis.text.x=element_text(size=8),
      axis.line = element_line(colour = 'black', size = 0.3)) +
ylab('Proportion of outbreak size') +
xlab('Median individual out-degree')


outbreak_dat %>%
  mutate(phase=cut(date,breaks=c(ymd('2020-01-01'),ymd('2020-10-15'),ymd('2021-01-01')),labels=c('Growth phase','Decline phase'))) %>%
  count(child,phase)

outbreak_dat %>%
  count(child,is.na(score)|score<0.95)

## Descriptive statistics
outbreak_dat %>% filter(!is_imputed)

outbreak_dat %>% count(q_type_binary)

get_dat_pairs(outbreak_dat) %>%
distinct(t_id)

outbreak_dat %>%
filter(!is.na(date_symptoms) & date_symptoms<date)

#hyperparameters
hyperparam_dat <- lapply(1:4,function(i){
                    readRDS(paste0(path_results,'outbreaker_res/',model_names[i],'.rds')) %>%
                    as_tibble() %>%
                    filter(step>=5000) %>%
                    mutate(chain=i) %>%
                    select(chain,step,eps,lambda,pi,mu)
                  }) %>% bind_rows()

pivot_longer(hyperparam_dat,cols=c('eps','lambda','pi','mu')) %>% 
group_by(name) %>%    
summarise(lower=quantile(value,0.025),
          mean=mean(value),
          upper=quantile(value,0.975)) %>%
write_xlsx(paste0(out_dir,'/hyperparameters.xlsx'))

s_i_flaxman <- tibble(s_i=1:100,p=get_serial_interval_distr_symptoms(N=100))

s_i_plot <- select(outbreaker_res_dat,iter,t_id,t_id_from,t_inf,kappa) %>% 
            inner_join(distinct(outbreaker_res_dat,iter,t_id,t_inf),by=c('iter','t_id_from'='t_id'),suffix=c('','_from')) %>% 
            mutate(s_i=t_inf-t_inf_from) %>% 
            filter(kappa==1) %>% 
            count(iter,s_i) %>%
            group_by(iter) %>% 
            mutate(p=n/sum(n)) %>%
            group_by(s_i) %>%
            summarise(mean_p=mean(p)) %>%
            ungroup() %>%
            left_join(s_i_flaxman,by='s_i') %>%
            mutate(p=if_else(is.na(p),0,p)) %>%
            ggplot() +
            geom_col(aes(s_i,mean_p),col='black',fill=my_colors(2)[2]) +
            geom_path(aes(s_i,p),col=my_colors(2)[1]) +
            xlab('Generation time') +
            ylab('Probability') +
            scale_x_continuous(limits=c(0.5,30),breaks = scales::pretty_breaks(6)) +
            scale_y_continuous(expand=c(0.001,0),breaks=scales::pretty_breaks(5),labels=scales::percent_format(2)) +
            theme_classic() +
            immunity_theme() +
            theme(axis.title.y=element_text(size=8),
                  axis.title.x=element_text(size=8),
                  axis.text.y=element_text(size=8),
                  axis.text.x=element_text(size=8),
                  axis.line = element_line(colour = 'black', size = 0.3))


inc_flaxman <- tibble(inc=1:100,p=get_incubation_distr_symptoms(N=100))
inc_plot <- filter(outbreaker_res_dat,!is.na(date_symptoms)) %>%
            mutate(t_inf_date=min(date)+t_inf) %>% 
            mutate(inc=as.numeric(date_symptoms-t_inf_date)) %>%
            count(iter,inc) %>%
            group_by(iter) %>% 
            mutate(p=n/sum(n)) %>%
            group_by(inc) %>%
            summarise(mean_p=mean(p)) %>%
            ungroup() %>%
            left_join(inc_flaxman,by='inc') %>%
            mutate(p=if_else(is.na(p),0,p)) %>%
            ggplot() +
            geom_col(aes(inc,mean_p),col='black',fill=my_colors(2)[2]) +
            geom_path(aes(inc,p),col=my_colors(2)[1]) +
            xlab('Incubation time') +
            ylab('Probability') +
            scale_x_continuous(limits=c(0.5,30),breaks = scales::pretty_breaks(6)) +
            scale_y_continuous(limits=c(0,0.21),expand=c(0.001,0),breaks=scales::pretty_breaks(5),labels=scales::percent_format(2)) +
            theme_classic() +
            immunity_theme() +
            theme(axis.title.y=element_text(size=8),
                  axis.title.x=element_text(size=8),
                  axis.text.y=element_text(size=8),
                  axis.text.x=element_text(size=8),
                  axis.line = element_line(colour = 'black', size = 0.3))

sup_figure_time_distr <- grid.arrange(inc_plot,s_i_plot,ncol=2)
ggsave(sup_figure_time_distr,filename = paste0(out_dir,'/sup_figure_time_distr.png'),width=180,height=60,units = 'mm',dpi = 320)
