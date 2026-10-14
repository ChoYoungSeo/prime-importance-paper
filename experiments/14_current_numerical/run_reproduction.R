#!/usr/bin/env Rscript
# The scientific functions are sourced unchanged. Only paths and run dispatch differ.
script<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE))
root<-normalizePath(file.path(dirname(script),'../..'));D<-file.path(root,'results/numerical_20261001')
args<-commandArgs(TRUE);stopifnot(length(args)==1L,args %in% c('smoke','all'))
for(pkg in c('digest','jsonlite','sensitivity','FNN','RANN','gtools','prime'))
 if(!requireNamespace(pkg,quietly=TRUE))stop('Missing dependency: ',pkg)
stopifnot(normalizePath(find.package('prime'))==normalizePath(file.path(root,'r-library/prime')),
 as.character(packageVersion('prime'))=='0.5.2.9101')
source(file.path(D,'reference_source/original_standardization_adapter.R'))
source(file.path(D,'reference_source/first_adapter.R'))
official<-new.env(parent=asNamespace('sensitivity'))
sys.source(file.path(D,'reference_source/sensitivity_1.27.1_pme.R'),envir=official)
first<-load_first_reference(file.path(root,'.verified/first_reference.R'))
runtime<-list(official=official,first=first)
lockEnvironment(official,bindings=TRUE);lockEnvironment(first,bindings=TRUE)
read<-function(p)read.csv(file.path(D,p),check.names=FALSE)
synthetic<-Filter(function(z)z$study=='synthetic',readRDS(file.path(D,'inputs/friedman_and_historical_concrete.rds')))
concrete<-readRDS(file.path(D,'inputs/concrete824.rds'))
stopifnot(length(synthetic)==270L,length(concrete)==101L)
if(args=='smoke'){
 synthetic<-Filter(function(z)z$task_id=='friedman_n100_rho0.9_r001',synthetic)
 concrete<-concrete['n0824_r01']
}
refF<-read('friedman/native_phi.csv');refL<-read('friedman/full_ll_allocations.csv')
refC<-rbind(read('concrete/paired_phi.csv'),read('concrete/full_data_reference.csv'))
falseF<-read('figures/figure4_native.csv');falseC<-read('figures/figure5_native.csv');fullC<-read('figures/figure5_references.csv')
inputsF<-read('friedman/input_manifest.csv');inputsC<-read('concrete/standardization_receipts.csv')
stamp<-paste(format(Sys.time(),'%Y%m%dT%H%M%S',tz='UTC'),Sys.getpid(),sep='-')
out<-file.path(root,'tmp/current_reproduction',stamp);dir.create(out,recursive=TRUE)
checks<-list();k<-0L
for(study in c('friedman','concrete'))for(z in if(study=='friedman')synthetic else concrete){
 id<-if(study=='friedman')z$task_id else z$id;s<-standardize_input(z$X)
 receipt<-if(study=='friedman')inputsF[inputsF$task_id==id,] else inputsC[inputsC$id==id,]
 expected_hash<-if(study=='friedman')receipt$standardized_X_hash else receipt$standardized_hash
 stopifnot(nrow(receipt)==1L,identical(s$standardized_hash,expected_hash))
 methods<-c('PRIME','PME-kNN','FIRST-SHAP',if(study=='friedman')'Full LL')
 # FALSE supplementary evidence is available only for the 30 displayed Friedman cases.
 if(study=='concrete'||(z$n==100&&z$rho==.9))methods<-c(methods,'PME-kNN FALSE')
 for(method in methods){
  seed<-if(study=='friedman')z$seed+switch(method,PRIME=100L,'PME-kNN'=200L,'FIRST-SHAP'=300L,'Full LL'=0L,'PME-kNN FALSE'=200L) else switch(method,PRIME=z$seed,'PME-kNN'=z$official_seed,'FIRST-SHAP'=z$first_seed,'PME-kNN FALSE'=z$official_seed)
  RNGkind('Mersenne-Twister','Inversion','Rejection');set.seed(seed)
  started<-proc.time()[['elapsed']]
  raw<-switch(method,
   PRIME=fit_prime_standardized(s$Z,z$y,runtime),
   'PME-kNN'=fit_pme_standardized(s$Z,z$y,runtime),
   'FIRST-SHAP'=fit_first_standardized(s$Z,z$y,runtime),
   'Full LL'={prepared<-prime:::.prime_prepare(s$Z,z$y,rescale=FALSE,twin_r=NULL,lambda=1e-6);fit<-prime:::.prime_estimate_subset_values(prepared,seq_len(ncol(s$Z)));list(effects=as.numeric(prime::prime_allocate(fit$raw)),fit=fit)},
   'PME-kNN FALSE'=official$pme_knn(model=function(X)z$y,X=s$Z,method='knn',n.knn=6,noise=FALSE,rescale=FALSE,marg=TRUE,tol=NULL,n.limit=2001,nboot=NULL,parl=NULL))
  elapsed<-proc.time()[['elapsed']]-started
  effects<-as.numeric(if(method %in% c('PME-kNN','PME-kNN FALSE'))raw$PME else raw$effects)
  if(method=='Full LL'){reference<-refL[refL$task_id==id&refL$rule=='PME',];reference<-reference[order(reference$variable_index),]}
  else if(method=='PME-kNN FALSE'){
   source<-if(study=='friedman')falseF else if(z$replication==0)fullC else falseC
   reference<-source[source$method=='PME-kNN (noise = FALSE)',]
   if('replication'%in%names(reference))reference<-reference[reference$replication==z$replication,]
   reference<-reference[order(reference$variable),]
  } else {
   reference<-if(study=='friedman')refF[refF$task_id==id&refF$method==method,] else refC[refC$id==id&refC$method==method,]
   reference<-reference[order(reference$variable),]
  }
  gap<-max(abs(effects-reference$estimate));tol<-1e-10
  passed<-length(effects)==ncol(z$X)&&nrow(reference)==length(effects)&&all(is.finite(effects))&&gap<=tol
  k<-k+1L;checks[[k]]<-data.frame(study=study,id=id,method=method,seed=seed,status=if(passed)'PASS' else 'FAIL',max_abs_difference=gap,tolerance=tol,elapsed_sec=elapsed)
  saveRDS(list(study=study,id=id,method=method,seed=seed,raw=raw,effects=effects,standardized_hash=s$standardized_hash),file.path(out,paste0(id,'.',gsub('[^A-Za-z0-9]','_',method),'.rds')))
  write.csv(do.call(rbind,checks),file.path(out,'parity.csv'),row.names=FALSE)
  cat(study,id,method,'gap=',format(gap,digits=16),'status=',if(passed)'PASS' else 'FAIL','\n')
  if(!passed)stop('Fresh output differs from saved native reference; inspect ',out)
 }
}
receipt<-list(status='PASS',scope=args,run_count=k,new_fits=k,method_source_bodies_unchanged=TRUE,reference='bundled native CSV estimates',comparison_tolerance=1e-10,grid_run_performed=args=='all',output=out,R=R.version.string,dependencies=as.list(setNames(vapply(c('prime','sensitivity','FNN','digest','RANN','gtools'),function(p)as.character(packageVersion(p)),''),c('prime','sensitivity','FNN','digest','RANN','gtools'))))
jsonlite::write_json(receipt,file.path(out,'RECEIPT.json'),pretty=TRUE,auto_unbox=TRUE)
cat('REPRODUCTION PASS: ',out,'\n')
