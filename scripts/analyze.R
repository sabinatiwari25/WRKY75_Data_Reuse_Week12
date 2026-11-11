source('scripts/setup_library.R')
suppressPackageStartupMessages({library(DESeq2); library(ggplot2); library(pheatmap); library(jsonlite)})
set.seed(75)
for (d in c('03_Results','04_Figures','06_Documentation')) dir.create(d,showWarnings=FALSE)
alpha <- 0.05
meta <- read.csv('02_Validation/samples_exp2.csv',check.names=FALSE)
raw <- read.csv('02_Validation/counts_exp2.csv',check.names=FALSE)
stopifnot(nrow(meta)==12, !anyDuplicated(meta$sample), !anyDuplicated(raw$gene_id))
counts <- as.matrix(raw[,meta$sample]); rownames(counts)<-raw$gene_id
stopifnot(!anyNA(counts),all(counts>=0),all(counts==round(counts)))
storage.mode(counts)<-'integer'
meta$genotype<-factor(meta$genotype,levels=c('Col-0','wrky75'))
meta$nitrogen<-factor(meta$nitrogen,levels=c('cntr','HN'))
rownames(meta)<-meta$sample
stopifnot(all(table(meta$genotype,meta$nitrogen)==3))
model<-model.matrix(~genotype*nitrogen,meta)
stopifnot(qr(model)$rank==4)
write.csv(model,'03_Results/design_matrix.csv')
keep <- rowSums(counts>=10)>=3
write.csv(data.frame(gene_id=rownames(counts),retained=keep,total_reads=rowSums(counts)), '03_Results/gene_filter.csv',row.names=FALSE)
dds<-DESeqDataSetFromMatrix(counts[keep,],meta,~genotype*nitrogen)
dds<-DESeq(dds,parallel=FALSE)
coefs<-resultsNames(dds)
writeLines(coefs,'03_Results/model_coefficients.txt')
interaction_name<-coefs[grepl('genotype.*nitrogen',coefs)]
genotype_name<-coefs[grepl('^genotype_',coefs)]
nitrogen_name<-coefs[grepl('^nitrogen_',coefs)]
stopifnot(length(interaction_name)==1,length(genotype_name)==1,length(nitrogen_name)==1)
res <- list(
 interaction=results(dds,name=interaction_name,alpha=alpha),
 nitrogen_Col0=results(dds,name=nitrogen_name,alpha=alpha),
 nitrogen_wrky75=results(dds,contrast=list(c(nitrogen_name,interaction_name)),alpha=alpha),
 genotype_control=results(dds,name=genotype_name,alpha=alpha),
 genotype_HN=results(dds,contrast=list(c(genotype_name,interaction_name)),alpha=alpha))
tables <- lapply(res,function(x){z<-as.data.frame(x);z$gene_id<-rownames(z);z<-z[,c('gene_id',setdiff(names(z),'gene_id'))];z[order(z$padj,na.last=TRUE),]})
for (nm in names(tables)) write.csv(tables[[nm]],paste0('03_Results/DESeq2_',nm,'.csv'),row.names=FALSE)
summary_table<-do.call(rbind,lapply(names(tables),function(nm){z<-tables[[nm]];s<-!is.na(z$padj)&z$padj<alpha;data.frame(contrast=nm,tested_pvalue=sum(!is.na(z$pvalue)),tested_adjusted=sum(!is.na(z$padj)),significant=sum(s),positive=sum(s&z$log2FoldChange>0),negative=sum(s&z$log2FoldChange<0))}))
write.csv(summary_table,'03_Results/contrast_summary.csv',row.names=FALSE)
normalized<-counts(dds,normalized=TRUE)
write.csv(data.frame(gene_id=rownames(normalized),normalized,check.names=FALSE),'03_Results/normalized_counts.csv',row.names=FALSE)
v<-varianceStabilizingTransformation(dds,blind=FALSE)
saveRDS(dds,'03_Results/dds.rds')
saveRDS(v,'03_Results/vst.rds')
qc<-data.frame(sample=meta$sample,genotype=meta$genotype,nitrogen=meta$nitrogen,
 library_reads=colSums(counts),detected_genes=colSums(counts>0),size_factor=sizeFactors(dds))
write.csv(qc,'03_Results/sample_qc.csv',row.names=FALSE)
palette<-c('Col-0'='#237A80','wrky75'='#B44748')
theme_set(theme_minimal(base_size=12)+theme(panel.grid.minor=element_blank(),plot.title=element_text(face='bold')))
save_plot<-function(name,p,w=8,h=5){ggsave(paste0('04_Figures/',name,'.png'),p,width=w,height=h,dpi=180,bg='white');ggsave(paste0('04_Figures/',name,'.pdf'),p,width=w,height=h,device=grDevices::cairo_pdf)}
p<-ggplot(qc,aes(sample,library_reads/1e6,fill=genotype))+geom_col()+scale_fill_manual(values=palette)+labs(title='Sequencing depth of the 12 selected samples',x=NULL,y='Assigned counts (millions)')+theme(axis.text.x=element_text(angle=60,hjust=1))
save_plot('library_sizes',p,9,5)
pc<-plotPCA(v,intgroup=c('genotype','nitrogen'),returnData=TRUE)
pv<-round(100*attr(pc,'percentVar'),1)
write.csv(pc,'03_Results/PCA_coordinates.csv',row.names=FALSE)
pc$replicate<-sub('.*_rep','',pc$name)
p<-ggplot(pc,aes(PC1,PC2,color=genotype,shape=nitrogen))+geom_point(size=4)+geom_text(aes(label=replicate),nudge_y=0.8,show.legend=FALSE,size=3)+scale_color_manual(values=palette)+labs(title='Experiment 2: principal component analysis',subtitle='Top 500 variable genes; labels indicate biological replicate',x=paste0('PC1 (',pv[1],'%)'),y=paste0('PC2 (',pv[2],'%)'))
save_plot('PCA',p)
png('04_Figures/sample_distance.png',width=1500,height=1200,res=180)
annotation<-meta[,c('genotype','nitrogen')]
pheatmap(as.matrix(dist(t(assay(v)))),annotation_col=annotation,annotation_row=annotation,main='Sample distances on variance-stabilized counts',fontsize=8)
dev.off()
z<-tables$interaction
z$significance<-ifelse(!is.na(z$padj)&z$padj<alpha,'FDR < 0.05','Other')
z$neglog<- -log10(pmax(z$padj,1e-300))
p<-ggplot(z[is.finite(z$neglog),],aes(log2FoldChange,neglog,color=significance))+geom_point(alpha=.55,size=1)+geom_hline(yintercept=-log10(alpha),linetype=2,color='grey40')+scale_color_manual(values=c('FDR < 0.05'='#B44748','Other'='#A7B1B5'))+labs(title='WRKY75 genotype by nitrogen interaction',subtitle='Positive: greater HN response in wrky75 than Col-0',x='Interaction log2 fold change (unshrunken)',y='-log10(BH adjusted P)')
save_plot('volcano_interaction',p)
png('04_Figures/dispersion.png',width=1400,height=1000,res=180); plotDispEsts(dds,main='DESeq2 dispersion estimates'); dev.off()

# Published lists: strip only a terminal transcript suffix; preserve source rows.
up<-read.csv('02_Validation/senescence_up_published.csv',check.names=FALSE)
down<-read.csv('02_Validation/senescence_down_published.csv',check.names=FALSE)
published<-rbind(transform(up,direction='up'),transform(down,direction='down'))
published$valid_agi<-toupper(as.character(published$valid_agi))=='TRUE'
write.csv(published,'03_Results/senescence_harmonized.csv',row.names=FALSE)
write.csv(published[!published$valid_agi,],'03_Results/senescence_unmapped.csv',row.names=FALSE)
sen_up<-unique(up$gene_id[toupper(as.character(up$valid_agi))=='TRUE'])
sen_down<-unique(down$gene_id[toupper(as.character(down$valid_agi))=='TRUE'])
sen_all<-union(sen_up,sen_down)
conflicts<-intersect(sen_up,sen_down)
write.csv(published[published$gene_id%in%conflicts,],'03_Results/senescence_direction_conflicts.csv',row.names=FALSE)
sig<-z$gene_id[!is.na(z$padj)&z$padj<alpha]
overlap<-intersect(sig,sen_all)
ol<-z[match(overlap,z$gene_id),]
ol$senescence_direction<-ifelse(ol$gene_id%in%conflicts,'both directions',ifelse(ol$gene_id%in%sen_up,'up in RNAi','down in RNAi'))
ol$senescence_source_rows<-vapply(ol$gene_id,function(g)sum(published$gene_id==g),integer(1))
# Keep transcript-level source effects as a range rather than choosing an arbitrary transcript.
ol$senescence_min_log2FC<-vapply(ol$gene_id,function(g)min(published$published_log2FC[published$gene_id==g]),numeric(1))
ol$senescence_max_log2FC<-vapply(ol$gene_id,function(g)max(published$published_log2FC[published$gene_id==g]),numeric(1))
write.csv(ol,'03_Results/shared_interaction_senescence.csv',row.names=FALSE)
descriptive<-do.call(rbind,lapply(names(tables),function(nm){q<-tables[[nm]];g<-q$gene_id[!is.na(q$padj)&q$padj<alpha];data.frame(contrast=nm,significant=length(g),shared_up=length(intersect(g,sen_up)),shared_down=length(intersect(g,sen_down)),shared_total=length(intersect(g,sen_all)))}))
write.csv(descriptive,'03_Results/descriptive_overlaps.csv',row.names=FALSE)
writeLines('No cross-study overlap enrichment test was performed: the complete tested-gene universe for GSE102302 is unavailable in the supplied significant-gene lists. A whole-genome or nitrogen-only universe would not establish comparable ascertainment in both studies.', '03_Results/overlap_universe_limitation.txt')
# Venn diagram is descriptive; circle areas are not proportional.
png('04_Figures/overlap.png',width=1500,height=1000,res=180)
par(mar=c(1,1,3,1));plot.new();plot.window(xlim=c(0,10),ylim=c(0,6),asp=1)
symbols(c(3.5,6.5),c(3,3),circles=c(2.2,2.2),inches=FALSE,add=TRUE,bg=c(adjustcolor('#237A80',.3),adjustcolor('#B44748',.3)),fg=c('#237A80','#B44748'))
text(c(2.8,5,7.2),3,c(length(setdiff(sig,sen_all)),length(overlap),length(setdiff(sen_all,sig))),cex=1.5)
text(2.7,.3,paste0('Nitrogen interaction\n',length(sig),' significant genes'));text(7.3,.3,paste0('Senescence published lists\n',length(sen_all),' unique mapped genes'))
title('Shared WRKY75-associated genes');mtext('Descriptive only; areas not proportional; no cross-study enrichment P-value',side=3,line=.2,cex=.75);dev.off()
heat_genes<-head(ol$gene_id,30)
heat_basis<-'Top shared genes ranked by interaction FDR'
if(length(heat_genes)<2){heat_genes<-head(sig,30);heat_basis<-'Top interaction genes ranked by interaction FDR (fewer than two shared genes)'}
if(length(heat_genes)>=2){
 png('04_Figures/expression_heatmap.png',width=1700,height=1600,res=180)
 pheatmap(assay(v)[heat_genes,],scale='row',annotation_col=annotation,cluster_cols=FALSE,main=heat_basis,fontsize=8,border_color=NA)
 dev.off()
 write.csv(data.frame(gene_id=heat_genes),'03_Results/heatmap_genes.csv',row.names=FALSE)
}
# GO BP over-representation within testable, annotated nitrogen genes.
go_status<-'not run'
go_summary<-list()
if(requireNamespace('org.At.tair.db',quietly=TRUE)&&requireNamespace('GO.db',quietly=TRUE)){
 suppressPackageStartupMessages({library(AnnotationDbi);library(org.At.tair.db);library(GO.db)})
 eligible<-z$gene_id[!is.na(z$padj)&grepl('^AT[1-5CM]G[0-9]{5}$',z$gene_id)]
 keys_valid<-intersect(eligible,keys(org.At.tair.db,keytype='TAIR'))
 mapping<-AnnotationDbi::select(org.At.tair.db,keys=keys_valid,keytype='TAIR',columns=c('GOALL','ONTOLOGYALL'))
 mapping<-unique(mapping[!is.na(mapping$GOALL)&mapping$ONTOLOGYALL=='BP',c('TAIR','GOALL')])
 universe<-unique(mapping$TAIR)
 query<-intersect(sig,universe)
 write.csv(mapping,'03_Results/GO_BP_mapping.csv',row.names=FALSE)
 write.csv(data.frame(gene_id=universe),'03_Results/GO_universe.csv',row.names=FALSE)
 terms<-split(mapping$TAIR,mapping$GOALL)
 terms<-terms[lengths(terms)>=10 & lengths(terms)<=500]
 if(length(query)>0&&length(query)<length(universe)&&length(terms)>0){
   go<-do.call(rbind,lapply(names(terms),function(id){g<-unique(terms[[id]]);k<-length(intersect(query,g));data.frame(GO=id,query_hits=k,term_universe=length(g),query_size=length(query),universe_size=length(universe),fold_enrichment=(k/length(query))/(length(g)/length(universe)),pvalue=phyper(k-1,length(g),length(universe)-length(g),length(query),lower.tail=FALSE),genes=paste(intersect(query,g),collapse=';'))}))
   go$padj<-p.adjust(go$pvalue,'BH')
   desc<-AnnotationDbi::select(GO.db,keys=go$GO,keytype='GOID',columns='TERM')
   go$description<-desc$TERM[match(go$GO,desc$GOID)]
   go<-go[order(go$padj,go$pvalue),]
   write.csv(go,'03_Results/GO_BP_interaction.csv',row.names=FALSE)
   significant_go<-go[go$padj<alpha,]
   go_status<-paste(nrow(significant_go),'BP terms at BH FDR < 0.05')
   go_summary<-list(universe=length(universe),query=length(query),terms_tested=nrow(go),significant=nrow(significant_go))
   if(nrow(significant_go)>0){
     top<-head(significant_go,15);top$description<-factor(top$description,levels=rev(top$description))
     p<-ggplot(top,aes(fold_enrichment,description,size=query_hits,color=-log10(padj)))+geom_point()+scale_color_gradient(low='#237A80',high='#B44748')+labs(title='GO Biological Process enrichment: interaction genes',subtitle='Annotated, testable nitrogen genes as background; GO terms of size 10-500',x='Fold enrichment',y=NULL,color='-log10(FDR)',size='Genes')
     save_plot('GO_BP',p,10,7)
   }
 } else {go_status<-'Not testable: empty or degenerate annotated query'}
}
cooks<-assays(dds)[['cooks']]
write.csv(data.frame(sample=colnames(dds),max_Cooks=apply(cooks,2,max,na.rm=TRUE)),'03_Results/Cooks_by_sample.csv',row.names=FALSE)
info<-list(input_genes=nrow(counts),retained_genes=sum(keep),excluded_low_count=sum(!keep),samples=ncol(counts),alpha=alpha,
 interaction_coefficient=interaction_name,coefficient_names=coefs,contrasts=summary_table,
 senescence=list(up_rows=nrow(up),down_rows=nrow(down),up_unique=length(sen_up),down_unique=length(sen_down),unique_total=length(sen_all),invalid_rows=sum(!published$valid_agi),direction_conflicts=length(intersect(sen_up,sen_down))),
 shared_interaction=length(overlap),overlap_test='not performed: common tested universe unavailable',
 PCA_variance_percent=pv,GO_status=go_status,GO=go_summary,heatmap_basis=heat_basis,
 library_size_range=range(qc$library_reads),size_factor_range=range(sizeFactors(dds)))
write_json(info,'03_Results/analysis_summary.json',pretty=TRUE,auto_unbox=TRUE,na='null')
writeLines(sub('[[:blank:]]+$','',capture.output(sessionInfo())), '06_Documentation/sessionInfo.txt')
packages<-sort(loadedNamespaces());write.csv(data.frame(package=packages,version=vapply(packages,function(p)as.character(getNamespaceVersion(p)),character(1))),'06_Documentation/package_versions.csv',row.names=FALSE)
writeLines(as.character(citation('DESeq2')),'06_Documentation/DESeq2_citation.txt')
print(summary_table);print(info)
