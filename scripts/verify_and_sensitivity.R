source('scripts/setup_library.R')
suppressPackageStartupMessages({library(DESeq2);library(jsonlite);library(AnnotationDbi);library(org.At.tair.db)})
dds<-readRDS('03_Results/dds.rds')
ids<-function(x)x$gene_id[!is.na(x$padj)&x$padj<0.05]
rd<-function(s)read.csv(paste0('03_Results/DESeq2_',s,'.csv'))
inter<-rd('interaction'); col_n<-rd('nitrogen_Col0'); mut_n<-rd('nitrogen_wrky75')
ctrl_g<-rd('genotype_control'); hn_g<-rd('genotype_HN')
lfc<-function(tab,rows)tab$log2FoldChange[match(rows,tab$gene_id)]
err_n<-max(abs(lfc(mut_n,inter$gene_id)-lfc(col_n,inter$gene_id)-inter$log2FoldChange),na.rm=TRUE)
err_g<-max(abs(lfc(hn_g,inter$gene_id)-lfc(ctrl_g,inter$gene_id)-inter$log2FoldChange),na.rm=TRUE)
stopifnot(err_n<1e-8,err_g<1e-8)
# DESeq2 BH adjustment applies to rows that pass its independent filter.
ok<-!is.na(inter$padj)
err_bh<-max(abs(p.adjust(inter$pvalue[ok],method='BH')-inter$padj[ok]))
stopifnot(err_bh<1e-10)
raw<-read.csv('02_Validation/counts_exp2.csv',check.names=FALSE)
m<-read.csv('02_Validation/samples_exp2.csv',check.names=FALSE)
stopifnot(identical(colnames(dds),m$sample),all(colSums(raw[,m$sample])==read.csv('03_Results/sample_qc.csv')$library_reads))
sen<-read.csv('03_Results/senescence_harmonized.csv')
sen_ids<-unique(sen$gene_id[sen$valid_agi]); shared<-read.csv('03_Results/shared_interaction_senescence.csv')
stopifnot(setequal(shared$gene_id,intersect(ids(inter),sen_ids)),!anyDuplicated(shared$gene_id))
go<-read.csv('03_Results/GO_BP_interaction.csv')
stopifnot(max(abs(go$padj-p.adjust(go$pvalue,'BH')))<1e-10)
stopifnot(max(abs(go$pvalue-phyper(go$query_hits-1,go$term_universe,go$universe_size-go$term_universe,go$query_size,lower.tail=FALSE)))<1e-10)

# Exploratory only: replicate index is visible in PCA but is not a documented block.
metadata<-as.data.frame(colData(dds));metadata$replicate_index<-factor(sub('.*_rep','',rownames(metadata)))
mm<-model.matrix(~replicate_index+genotype*nitrogen,metadata)
stopifnot(qr(mm)$rank==6)
blocked<-DESeqDataSetFromMatrix(counts(dds),metadata,~replicate_index+genotype*nitrogen)
blocked<-DESeq(blocked,parallel=FALSE)
nm<-resultsNames(blocked)[grepl('genotype.*nitrogen',resultsNames(blocked))]
br<-as.data.frame(results(blocked,name=nm,alpha=.05));br$gene_id<-rownames(br)
br<-br[order(br$padj,na.last=TRUE),c('gene_id','baseMean','log2FoldChange','lfcSE','stat','pvalue','padj')]
write.csv(br,'03_Results/sensitivity_replicate_index.csv',row.names=FALSE)
both<-intersect(ids(inter),ids(br)); both_shared<-intersect(both,sen_ids)
sgn<-sum(sign(lfc(inter,both))==sign(lfc(br,both)))
pc<-read.csv('03_Results/PCA_coordinates.csv');pc$replicate_index<-as.integer(sub('.*_rep','',pc$name))
group_order<-tapply(seq_len(nrow(pc)),pc$group,function(i)all(diff(pc$PC2[i[order(pc$replicate_index[i])]])>0)||all(diff(pc$PC2[i[order(pc$replicate_index[i])]])<0))
comparable<-is.finite(inter$log2FoldChange)&is.finite(lfc(br,inter$gene_id))
out<-list(contrast_max_abs_error_nitrogen=err_n,contrast_max_abs_error_genotype=err_g,BH_max_abs_error=err_bh,
    PC2_monotonic_groups=sum(group_order),PC2_groups=length(group_order),primary_significant=length(ids(inter)),
    sensitivity_significant=length(ids(br)),significant_in_both=length(both),sign_agreement_in_both=sgn,
    shared_senescence_in_both=length(both_shared),sensitivity_shared=length(intersect(ids(br),sen_ids)),
    interaction_log2FC_pearson=cor(inter$log2FoldChange[comparable],lfc(br,inter$gene_id)[comparable]))
write_json(out,'03_Results/verification.json',pretty=TRUE,auto_unbox=TRUE)
write.csv(data.frame(gene_id=both_shared),'03_Results/sensitivity_shared_robust.csv',row.names=FALSE)
# Annotation for a compact, auditable set of example loci in the report.
candidate<-unique(c(head(shared$gene_id,10),'AT5G13080','AT1G74710','AT4G35090','AT2G14610'))
annotation<-AnnotationDbi::select(org.At.tair.db,keys=candidate,keytype='TAIR',columns=c('SYMBOL','GENENAME'))
write.csv(annotation,'03_Results/candidate_annotations.csv',row.names=FALSE)
candidate_data<-merge(data.frame(gene_id=candidate),inter,by='gene_id',all.x=TRUE,sort=FALSE)
candidate_data$nitrogen_Col0_log2FC<-lfc(col_n,candidate_data$gene_id)
candidate_data$nitrogen_wrky75_log2FC<-lfc(mut_n,candidate_data$gene_id)
candidate_data$leaf_list<-ifelse(candidate_data$gene_id%in%sen_ids,'listed','not listed')
write.csv(candidate_data,'03_Results/candidate_statistics.csv',row.names=FALSE)
print(out)
