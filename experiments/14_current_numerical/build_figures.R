#!/usr/bin/env Rscript
# Rebuild numerical displays from saved estimates; never refit or modify native data.
suppressPackageStartupMessages(library(ggplot2))
script <- sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value=TRUE))
root <- normalizePath(file.path(dirname(script),'../..'))
data <- file.path(root,'results/numerical_20261001')
out <- file.path(root,'tmp/numerical_20261001'); dir.create(out,recursive=TRUE,showWarnings=FALSE)
status <- system2('python3',shQuote(file.path(root,'experiments/14_current_numerical/validate.py')),stdout=FALSE)
stopifnot(status==0L)
read <- function(p) read.csv(file.path(data,p),check.names=FALSE)
methods <- c('PRIME','PME-kNN (noise = TRUE)','FIRST-SHAP')
labels <- setNames(c('PRIME','PME-NA','FIRST-SHAP'),methods)
pal <- setNames(c('#24557A','#B37038','#6D697C'),methods)
t <- read('timing/summary.csv'); tm <- c(PRIME='PRIME',PME_kNN_TRUE=methods[2],FIRST_SHAP='FIRST-SHAP')
t <- t[t$method %in% names(tm),];t$method <- factor(unname(tm[t$method]),levels=methods)
for (num in 4:5) {
 d<-read(sprintf('figures/figure%d_native.csv',num));d<-d[d$method %in% methods,]
 d$method<-factor(d$method,levels=methods)
 p<-if(num==4)10L else 8L
 vnames<-if(num==4)paste0('X',1:10) else c('Cement','Slag','Fly ash','Water/binder','Superplasticizer','Coarse aggregate','Fine aggregate','Age')
 d$x<-factor(d$variable,levels=seq_len(p),labels=vnames)
 g<-ggplot(d,aes(x,estimate,fill=method))+geom_hline(yintercept=0,colour='grey82',linewidth=.25)+
  geom_boxplot(width=.64,linewidth=.3,outlier.shape=21,outlier.fill='white',outlier.size=.95)+
  facet_wrap(~method,nrow=1,labeller=as_labeller(labels))+scale_fill_manual(values=pal)+
  labs(x=NULL,y='Importance measure')+theme_classic(base_size=12)+
  theme(legend.position='none',strip.background=element_blank(),strip.text=element_text(face='bold'),
        axis.text.x=element_text(angle=if(num==4)0 else 43,hjust=if(num==4).5 else 1))
 if(num==4){
  refs<-read('figures/figure4_references.csv');refs<-refs[refs$method %in% methods,]
  refs$method<-factor(refs$method,levels=methods);refs$x<-factor(refs$variable,levels=seq_len(p),labels=vnames)
  g<-g+geom_point(data=refs,aes(y=estimate),shape=23,fill='white',colour='#4F7D64',size=2)
 }
 full<-g; detail<-g+coord_cartesian(ylim=c(-.30,1.05),expand=FALSE)
 stopifnot(identical(ggplot_build(full)$data[[2]],ggplot_build(detail)$data[[2]]))
 if(num==4){
  t$x<-5.5;t$y<-.99;t$label<-sprintf('%.3f s',t$median_sec)
  detail<-detail+geom_text(data=t,aes(x=x,y=y,label=label),inherit.aes=FALSE,colour='#C2312A')
 }
 for(ext in c('pdf','png')){
  ggsave(file.path(out,sprintf('figure%d_detail.%s',num,ext)),detail,width=12.4,height=if(num==4)3.8 else 4.8,dpi=180)
  ggsave(file.path(out,sprintf('figure%d_full.%s',num,ext)),full,width=12.4,height=if(num==4)3.8 else 4.8,dpi=180)
 }
}
cat('PASS: saved-output figure rebuild; full and detail box statistics identical. Author manuscript PDFs remain unchanged.\n')
