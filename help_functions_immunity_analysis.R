read_immunity_files <- function(subtrees,replay,init_size,efficacy,N,max_chains,max_sim,path){
  lapply(1:max_chains,function(ch){
    lapply(0:max_sim,function(i){
      fread(paste0(path,'immunity_res_',if_else(replay,'','no'),'replay_',if_else(subtrees,'','no'),'subtrees','_N',N,'_init',init_size,'_eff',efficacy,'_',ch,'_',i,'.tsv')) %>%
        mutate(sim=161*(ch-1) + i*N + sim)
    }) %>% bind_rows()
  }) %>% bind_rows()
}

plot_immunity <- function(immunity_dat){
  line_thickness <- 0.6
  point_size <- 1.5
  plot_dat <- immunity_dat %>%
    group_by(Strategy,p) %>%
    summarise(lower=quantile(mean,0.025),
              upper=quantile(mean,0.975),
              mean=mean(mean)) %>%
    ungroup()%>%
    #filter(p!=0.3) %>%
    mutate(Strategy=factor(Strategy,levels=c('Age down','Age up','Uniform'),labels=c('Age descending','Age ascending','Uniform at random')))
  ggplot(plot_dat) +
    geom_point(aes(x=p,y=mean,col=Strategy),size=point_size) +
    geom_line(aes(x=p,y=mean,col=Strategy),size=line_thickness) +
    geom_ribbon(aes(x=p,ymin=lower,ymax=upper,fill=Strategy),alpha=0.2) +
    geom_abline(intercept=100,slope=0,col='red',linetype='dashed',size=line_thickness) +
    immunity_theme() +
    ylab('Estimated size of third wave') +
    xlab('Percentage vaccinated')+
    scale_y_continuous(labels=c(0,100,500,1000,1500,2000,2500),breaks=c(0,100,500,1000,1500,2000,2500))+
    scale_x_continuous(labels=scales::percent_format(2),breaks=scales::pretty_breaks(6))
}

plot_immunity_prior <- function(immunity_dat){
  line_thickness <- 0.6
  point_size <- 1.5
  plot_dat <- immunity_dat %>%
    group_by(Strategy,p) %>%
    summarise(lower=quantile(mean,0.025),
              upper=quantile(mean,0.975),
              mean=mean(mean)) %>%
    ungroup()%>%
    filter(p!=0.3) %>%
    mutate(Strategy=factor(Strategy,levels=c('Replay','Age down','Age up','Uniform'),labels=c('Replay','Age descending','Age ascending','Uniform at random')))
  ggplot(plot_dat) +
    geom_point(data=filter(plot_dat,Strategy=='Replay'),aes(x=p,y=mean),size=point_size) +
    geom_line(data=filter(plot_dat,Strategy=='Replay'),aes(x=p,y=mean),size=line_thickness) +
    geom_ribbon(data=filter(plot_dat,Strategy=='Replay'),aes(x=p,ymin=lower,ymax=upper),alpha=0.2)+
    geom_point(data=filter(plot_dat,Strategy!='Replay'),aes(x=p,y=mean,col=Strategy),size=point_size) +
    geom_line(data=filter(plot_dat,Strategy!='Replay'),aes(x=p,y=mean,col=Strategy),size=line_thickness) +
    geom_ribbon(data=filter(plot_dat,Strategy!='Replay'),aes(x=p,ymin=lower,ymax=upper,fill=Strategy),alpha=0.2)+
    geom_abline(intercept=100,slope=0,col='red',linetype='dashed',size=line_thickness) +
    geom_vline(xintercept=0.295,linetype='dashed',size=line_thickness) +
    annotate('text',x=c(0.32,0.32),y=c(2000,1800),label=c('Vaccinations ','as of April 28th'),size=4,hjust=0) + 
    immunity_theme() +
    ylab('Estimated size of third wave') +
    xlab('Percentage vaccinated')+
    scale_y_continuous(labels=c(0,100,500,1000,1500,2000,2500),breaks=c(0,100,500,1000,1500,2000,2500))+
    scale_x_continuous(labels=scales::percent_format(2),breaks=scales::pretty_breaks(6))
}

plot_immunity_subtrees <- function(immunity_dat){
  line_thickness <- 0.6
  point_size <- 1.5
  plot_dat <- immunity_dat %>%
    group_by(Strategy,p) %>%
    summarise(lower=quantile(mean,0.025),
              upper=quantile(mean,0.975),
              mean=mean(mean)) %>%
    ungroup()%>%
    mutate(Strategy=factor(Strategy,levels=c('Age down','Age up','Uniform'),labels=c('Age descending','Age ascending','Uniform at random')))
  ggplot(plot_dat) +
    geom_point(aes(x=p,y=mean,col=Strategy),size=point_size) +
    geom_line(aes(x=p,y=mean,col=Strategy),size=line_thickness) +
    geom_ribbon(aes(x=p,ymin=lower,ymax=upper,fill=Strategy),alpha=0.2) +
    immunity_theme() +
    ylab('Percentage of original outbreak size') +
    xlab('Percentage vaccinated')+
    scale_y_continuous(limits=c(0,1.05),labels=scales::percent_format(2),breaks=scales::pretty_breaks(6),expand=c(0.01,0.01))+
    scale_x_continuous(labels=scales::percent_format(2),breaks=scales::pretty_breaks(6),expand=c(0.01,0.01))
}

plot_immunity_subtrees_prior <- function(immunity_dat){
  line_thickness <- 0.6
  point_size <- 1.5
  plot_dat <- immunity_dat %>%
    group_by(Strategy,p) %>%
    summarise(lower=quantile(mean,0.025),
              upper=quantile(mean,0.975),
              mean=mean(mean)) %>%
    ungroup()%>%
    filter(p!=0.3) %>%
    mutate(Strategy=factor(Strategy,levels=c('Replay','Age down','Age up','Uniform'),labels=c('Replay','Age descending','Age ascending','Uniform at random')))
  ggplot(plot_dat) +
    geom_point(data=filter(plot_dat,Strategy=='Replay'),aes(x=p,y=mean),size=point_size) +
    geom_line(data=filter(plot_dat,Strategy=='Replay'),aes(x=p,y=mean),size=line_thickness) +
    geom_ribbon(data=filter(plot_dat,Strategy=='Replay'),aes(x=p,ymin=lower,ymax=upper),alpha=0.2)+
    geom_point(data=filter(plot_dat,Strategy!='Replay'),aes(x=p,y=mean,col=Strategy),size=point_size) +
    geom_line(data=filter(plot_dat,Strategy!='Replay'),aes(x=p,y=mean,col=Strategy),size=line_thickness) +
    geom_ribbon(data=filter(plot_dat,Strategy!='Replay'),aes(x=p,ymin=lower,ymax=upper,fill=Strategy),alpha=0.2)+
    geom_vline(xintercept=0.295,linetype='dashed',size=line_thickness) +
    annotate('text',x=c(0.32,0.32),y=c(0.77,0.69),label=c('Vaccinations ','as of April 28th'),size=4,hjust=0) +
    immunity_theme() +
    ylab('Percentage of original outbreak size') +
    xlab('Percentage vaccinated')+
    scale_y_continuous(limits=c(0,1.05),labels=scales::percent_format(2),breaks=scales::pretty_breaks(6),expand=c(0.01,0.01))+
    scale_x_continuous(labels=scales::percent_format(2),breaks=scales::pretty_breaks(6),expand=c(0.01,0.01))
}

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
