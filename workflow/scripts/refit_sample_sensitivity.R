#!/usr/bin/env Rscript
# Refit the recorded expression models on two fixed sample sets.
# Counts, gene universe and functional-set membership are held constant.
options(stringsAsFactors=FALSE, warn=1)
suppressPackageStartupMessages({library(edgeR); library(limma); library(jsonlite)})
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=3L) stop('Usage: refit_sample_sensitivity.R INPUT_DIR OUTPUT_DIR SETTINGS_JSON')
in_dir <- args[1]; out_dir <- args[2]; cfg <- fromJSON(args[3])
read_tsv <- function(p, ...) read.delim(p, sep='\t', check.names=FALSE, quote='', comment.char='', ...)
write_tsv <- function(x,p) {
  dir.create(dirname(p), recursive=TRUE, showWarnings=FALSE)
  con <- if(grepl('\\.gz$',p)) gzfile(p,'wt') else p
  if(inherits(con,'connection')) on.exit(close(con))
  write.table(x,con,sep='\t',quote=FALSE,row.names=FALSE,col.names=TRUE,na='NA')
}
counts <- as.matrix(read_tsv(file.path(in_dir,'counts.tsv'),row.names=1))
if(any(!is.finite(counts)) || any(counts<0) || any(counts!=floor(counts))) stop('Invalid counts')
storage.mode(counts)<-'integer'
primary <- read_tsv(file.path(in_dir,'samples.tsv'))
genes <- read_tsv(file.path(in_dir,'genes.tsv'))$gene_id
if(anyDuplicated(genes) || !all(genes %in% rownames(counts))) stop('Invalid fixed gene universe')
counts <- counts[genes,,drop=FALSE]
sets <- read_tsv(file.path(in_dir,'sets.tsv'))
sets <- sets[sets$module!='hypothetical_uncharacterized',,drop=FALSE]
sets <- lapply(split(sets$gene_id,sets$module),unique)
sets <- lapply(sets,function(ids) intersect(ids,genes))
sets <- sets[lengths(sets)>=5 & lengths(sets)<=500]
if(!(cfg$focus_module %in% names(sets))) stop('Focus set missing')
if(!setequal(sets[[cfg$focus_module]],cfg$focus_genes)) stop('Focus membership changed')
levels_g <- c('wt_dry','bmr6_dry','bmr12_dry','wt_wet','bmr6_wet','bmr12_wet')
contrast_matrix <- function(cols) makeContrasts(
  bmr6_vs_wt_dry=bmr6_dry-wt_dry,
  bmr12_vs_wt_dry=bmr12_dry-wt_dry,
  bmr6_vs_wt_wet=bmr6_wet-wt_wet,
  bmr12_vs_wt_wet=bmr12_wet-wt_wet,
  wet_vs_dry_wt=wt_wet-wt_dry,
  wet_vs_dry_bmr6=bmr6_wet-bmr6_dry,
  wet_vs_dry_bmr12=bmr12_wet-bmr12_dry,
  bmr6_water_interaction=(bmr6_wet-bmr6_dry)-(wt_wet-wt_dry),
  bmr12_water_interaction=(bmr12_wet-bmr12_dry)-(wt_wet-wt_dry),
  bmr12_vs_bmr6_dry=bmr12_dry-bmr6_dry,
  bmr12_vs_bmr6_wet=bmr12_wet-bmr6_wet,levels=cols)
run_model <- function(samples,name,scale) {
  if(anyDuplicated(samples$sample_id)) stop('Duplicate samples')
  samples$group<-factor(samples$group,levels=levels_g)
  rownames(samples)<-samples$sample_id
  if(anyNA(samples$group) || any(table(samples$group)<2)) stop('Invalid or under-replicated cells')
  x <- counts[,samples$sample_id,drop=FALSE]
  if(scale=='state') {
    y<-DGEList(counts=x,samples=samples)
    y<-calcNormFactors(y,method='TMMwsp')
  } else {
    lib<-as.numeric(samples$all_primary)
    if(any(!is.finite(lib)) || any(lib<=0)) stop('Invalid tissue offsets')
    y<-DGEList(counts=x,samples=samples,lib.size=lib)
    y$samples$norm.factors<-1
  }
  design<-model.matrix(~0+group,data=samples)
  colnames(design)<-sub('^group','',colnames(design)); rownames(design)<-samples$sample_id
  if(qr(design)$rank!=ncol(design)) stop('Rank-deficient design')
  y<-estimateDisp(y,design,robust=TRUE)
  fit<-glmQLFit(y,design,robust=TRUE)
  cm<-contrast_matrix(colnames(design)); dest<-file.path(out_dir,name)
  logcpm<-cpm(y,log=TRUE,prior.count=2)
  write_tsv(data.frame(gene_id=rownames(logcpm),logcpm,check.names=FALSE),file.path(dest,'logCPM.tsv.gz'))
  write_tsv(data.frame(sample_id=colnames(y),library_size=y$samples$lib.size,
       normalization_factor=y$samples$norm.factors,
       effective_library_size=y$samples$lib.size*y$samples$norm.factors),file.path(dest,'normalization_factors.tsv'))
  write_tsv(samples,file.path(dest,'samples.tsv'))
  write_tsv(data.frame(sample_id=rownames(design),design,check.names=FALSE),file.path(dest,'design.tsv'))
  indices<-ids2indices(sets,rownames(y),remove.empty=TRUE)
  v<-voom(y,design=design,plot=FALSE,normalize.method='none')
  cams<-list(); des<-list()
  for(ct in colnames(cm)) {
    test<-glmQLFTest(fit,contrast=cm[,ct]); tab<-topTags(test,n=Inf,sort.by='none')$table
    tab$gene_id<-rownames(tab); rownames(tab)<-NULL
    write_tsv(tab[,c('gene_id',setdiff(names(tab),'gene_id'))],file.path(dest,'de',paste0(ct,'.tsv.gz')))
    des[[ct]]<-data.frame(contrast=ct,genes=nrow(tab),fdr05=sum(tab$FDR<0.05),
             fdr05_lfc1=sum(tab$FDR<0.05 & abs(tab$logFC)>=1))
    # Preserve the archived CAMERA calculation: voom E matrix, not voom precision weights.
    cam<-camera(v$E,index=indices,design=design,contrast=cm[,ct],
                inter.gene.cor=0.01,sort=FALSE)
    cam$module<-rownames(cam); cam$contrast<-ct; rownames(cam)<-NULL; cams[[ct]]<-cam
  }
  camera_all<-do.call(rbind,cams); rownames(camera_all)<-NULL
  camera_all$FDR_all_set_contrast_pairs<-p.adjust(camera_all$PValue,method='BH')
  write_tsv(camera_all,file.path(dest,'camera.tsv'))
  write_tsv(do.call(rbind,des),file.path(dest,'contrast_summary.tsv'))
  info<-list(model=name,n_samples=ncol(y),n_genes=nrow(y),n_sets=length(sets),
    normalization=if(scale=='state') 'TMMwsp' else 'all_primary_post_filter_offset',
    prior_count=2,robust_dispersion=TRUE,robust_ql=TRUE,camera_inter_gene_cor=0.01,
    camera_uses_voom_weights=FALSE,camera_FDR_family='sets_within_each_contrast',
    common_dispersion=unname(y$common.dispersion))
  write_json(info,file.path(dest,'model.json'),pretty=TRUE,auto_unbox=TRUE)
  message('MODEL_COMPLETE=',name)
}
if(nrow(primary)!=cfg$expected_primary_n) stop('Unexpected primary sample count')
if(!all(cfg$exclude_samples %in% primary$sample_id)) stop('Exclusion sample absent')
reduced<-primary[!primary$sample_id %in% cfg$exclude_samples,,drop=FALSE]
if(nrow(reduced)!=cfg$expected_reduced_n) stop('Unexpected reduced sample count')
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
write_tsv(data.frame(package=c('R','edgeR','limma','jsonlite','statmod'),
  version=c(R.version.string,as.character(packageVersion('edgeR')),as.character(packageVersion('limma')),
            as.character(packageVersion('jsonlite')),as.character(packageVersion('statmod')))),file.path(out_dir,'versions.tsv'))
writeLines(capture.output(sessionInfo()),file.path(out_dir,'sessionInfo.txt'))
for(scale in c('state','tissue')) {
  run_model(primary,paste0('baseline22_',scale),scale)
  run_model(reduced,paste0('exclude21_',scale),scale)
}
writeLines('complete',file.path(out_dir,'REFIT_COMPLETE.txt'))
