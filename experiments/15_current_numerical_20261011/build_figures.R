#!/usr/bin/env Rscript
# Portable version of the approved 2026-10-11 R plotting layers. No fitting.
args <- commandArgs(trailingOnly=TRUE)
arg <- function(flag,default=NULL) {i<-match(flag,args);if(is.na(i))return(default);if(i==length(args))stop('Missing value for ',flag);args[i+1L]}
dataset <- arg('--dataset');stopifnot(dataset %in% c('friedman','concrete'))
script <- sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE))
root <- normalizePath(file.path(dirname(script),'../..'))
data_root <- file.path(root,'results/numerical_20261011')
out <- arg('--output-dir',file.path(root,'build/current_numerical_20261011'))
dir.create(out,recursive=TRUE,showWarnings=FALSE);out<-normalizePath(out)
stopifnot(!startsWith(paste0(out,'/'),paste0(file.path(root,'manuscript'),'/')),
          !startsWith(paste0(out,'/'),paste0(file.path(root,'results'),'/')))
.libPaths(c(file.path(root,'r-library'),.libPaths()))
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(jsonlite)})
# Grob measurement needs a device; avoid R's implicit Rplots.pdf in the checkout.
grDevices::pdf(file=NULL)
num<-if(dataset=='friedman')4L else 5L;p<-if(num==4L)10L else 8L
R_expected<-if(num==4L)30L else 100L
vnames<-if(num==4L)paste0('X',1:10) else c('Cement','Slag','Fly ash','Water/binder','Superplasticizer','Coarse aggregate','Fine aggregate','Age')
old_methods<-c('PRIME','PME-kNN (noise = TRUE)','FIRST-SHAP')
old_labels<-setNames(c('PRIME','PME-NA','FIRST-SHAP'),old_methods)
old_palette<-setNames(c('#24557A','#B37038','#6D697C'),old_methods)
native_method<-function(x) ifelse(x=='PME-NA','PME-kNN (noise = TRUE)',x)
all_values<-read.csv(file.path(data_root,'results/plot_values.csv'),check.names=FALSE)
d<-all_values[all_values$dataset==dataset & all_values$method!='SHAP',]
d$method<-factor(native_method(d$method),levels=old_methods)
d$x<-factor(d$variable,levels=seq_len(p),labels=vnames)
refs<-read.csv(file.path(data_root,'results/friedman_references.csv'),check.names=FALSE)
refs$method<-factor(native_method(refs$method),levels=old_methods)
refs$x<-factor(refs$variable,levels=seq_len(p),labels=vnames)
timing<-read.csv(file.path(data_root,'results/timing_summary.csv'),check.names=FALSE)
timing<-timing[timing$dataset==dataset,]
labels<-setNames(timing$display,timing$method)
z<-all_values[all_values$dataset==dataset & all_values$method=='SHAP',]
z$x<-factor(z$variable,levels=seq_len(p),labels=vnames)
stopifnot(nrow(d)==3L*p*R_expected,nrow(z)==p*R_expected)
original_plot<-function(dat){
 ggplot(dat,aes(x,estimate,fill=method))+geom_hline(yintercept=0,colour='grey82',linewidth=.25)+
 geom_boxplot(width=.64,linewidth=.3,outlier.shape=21,outlier.fill='white',outlier.size=.95)+
 facet_wrap(~method,nrow=1,labeller=as_labeller(old_labels))+scale_fill_manual(values=old_palette)+guides(fill='none')+
 labs(x=NULL,y='Importance measure')+theme_classic(base_size=12)+
 theme(legend.position='none',strip.background=element_blank(),strip.text=element_text(face='bold'),
 axis.text.x=element_text(angle=if(num==4L)0 else 43,hjust=if(num==4L).5 else 1))
}
canonical<-function(layer){cols<-c('x','ymin','lower','middle','upper','ymax','outliers');x<-layer[order(as.numeric(layer$x)),cols,drop=FALSE];rownames(x)<-NULL;x$x<-as.numeric(x$x);x}
approved<-fromJSON(file.path(data_root,'provenance/approved_figure_statistics.json'),simplifyVector=FALSE)[[dataset]]
as_rows<-function(x)fromJSON(toJSON(x,dataframe='rows',digits=17),simplifyVector=FALSE)
plots<-list();stats<-list()
for(i in seq_along(old_methods)){
 method<-old_methods[i];zz<-d[d$method==method,];g<-original_plot(zz)
 if(num==4L)g<-g+geom_point(data=refs[refs$method==method,],aes(y=estimate),shape=23,fill='white',colour='#4F7D64',size=2)
 y_view<-if(i==3L)c(-.05,1.05) else c(-.30,1.05)
 g<-g+coord_cartesian(ylim=y_view,expand=if(i==3L)c(top=FALSE,right=TRUE,bottom=FALSE,left=TRUE) else FALSE)
 label<-labels[[unname(old_labels[method])]]
 if(num==4L){labeldata<-data.frame(x=5.5,y=.99,label=label);g<-g+geom_text(data=labeldata,aes(x=x,y=y,label=label),inherit.aes=FALSE,colour='#C2312A')}
 if(num==5L)g<-g+annotate('text',x=(p+2)/2,y=Inf,label=label,vjust=1.4,colour='#C2312A')
 now<-canonical(ggplot_build(g)$data[[2]])
 stopifnot(isTRUE(all.equal(as_rows(now),approved$baseline_preservation[[method]]$box_statistics,check.attributes=FALSE,tolerance=1e-14)))
 stats[[unname(old_labels[method])]]<-now;plots[[i]]<-g
}
g<-ggplot(z,aes(x,estimate,group=x))+geom_hline(yintercept=0,colour='grey82',linewidth=.25)+
 geom_boxplot(fill='#737A82',width=.64,linewidth=.3,outlier.shape=21,outlier.fill='white',outlier.size=.95)+
 facet_wrap(~method,nrow=1,labeller=as_labeller(c(SHAP='Mean absolute SHAP')))+
 labs(x=NULL,y='Mean absolute SHAP values')+theme_classic(base_size=12)+
 theme(strip.background=element_blank(),strip.text=element_text(face='bold'),axis.text.x=element_text(angle=if(num==4L)0 else 43,hjust=if(num==4L).5 else 1),legend.position='none')
shap_before<-ggplot_build(g);shap_stats<-canonical(shap_before$data[[2]])
shap_original_range<-shap_before$layout$panel_params[[1]]$y.range
g<-g+coord_cartesian(ylim=shap_original_range,expand=c(top=FALSE,right=TRUE,bottom=FALSE,left=TRUE))
g<-g+annotate('text',x=if(num==4L)(p+1)/2 else (p+2)/2,y=Inf,label=labels[['SHAP']],vjust=1.4,colour='#C2312A')
stopifnot(identical(shap_stats,canonical(ggplot_build(g)$data[[2]])),
 isTRUE(all.equal(as_rows(shap_stats),approved$new_panel$box_statistics,check.attributes=FALSE,tolerance=1e-14)))
stats[['SHAP']]<-shap_stats;plots[[4]]<-g
figure<-wrap_plots(plots,ncol=2,nrow=2)+plot_layout(heights=c(1,1))&theme(legend.position='none')
first_range<-ggplot_build(plots[[3]])$layout$panel_params[[1]]$y.range
shap_range<-ggplot_build(plots[[4]])$layout$panel_params[[1]]$y.range
first_zero_npc <- -first_range[1]/diff(first_range)
shap_zero_npc <- -shap_range[1]/diff(shap_range)
stopifnot(first_range[1]==-.05,first_range[2]==1.05,identical(shap_range,shap_original_range),abs(first_zero_npc-shap_zero_npc)<1e-14)
first_xrange<-ggplot_build(plots[[3]])$layout$panel_params[[1]]$x.range
shap_xrange<-ggplot_build(plots[[4]])$layout$panel_params[[1]]$x.range
stopifnot(isTRUE(all.equal(first_xrange,c(.4,p+.6),tolerance=1e-14)),identical(first_xrange,shap_xrange))
figure_gtable<-patchwork::patchworkGrob(figure)
panel_layout<-figure_gtable$layout[grepl('^panel',figure_gtable$layout$name)&figure_gtable$layout$name!='panel-area',c('name','t','b','l','r')]
low<-panel_layout[grepl('-[34]$',panel_layout$name),]
stopifnot(nrow(low)==2L,length(unique(low$t))==1L,length(unique(low$b))==1L)
receipt<-list(status='PASS',dataset=dataset,replications=R_expected,all_four_box_statistics_match_approved=TRUE,
 runtime_labels=as.list(labels),zero_height_fraction=unname(first_zero_npc),lower_panel_same_gtable_rows=TRUE,
 first_SHAP_y_range=first_range,mean_absolute_SHAP_y_range=shap_range,lower_x_range=first_xrange,
 box_statistics=stats,packages=setNames(lapply(c('ggplot2','patchwork','jsonlite'),function(x)as.character(packageVersion(x))),c('ggplot2','patchwork','jsonlite')),
 no_fitting=TRUE,no_manuscript_files_written=TRUE)
if(!'--check'%in%args){
 stem<-file.path(out,if(num==4L)'PRIME_Friedman' else 'PRIME_Concrete')
 ggsave(paste0(stem,'.pdf'),figure_gtable,width=8.4,height=if(num==4L)7.6 else 9.6,bg='white')
 ggsave(paste0(stem,'.png'),figure_gtable,width=8.4,height=if(num==4L)7.6 else 9.6,dpi=180,bg='white')
}
write_json(receipt,file.path(out,paste0(dataset,'_build_validation.json')),pretty=TRUE,auto_unbox=TRUE,digits=17)
cat('PASS:',dataset,'approved four-panel statistics, labels and zero alignment.\n')
