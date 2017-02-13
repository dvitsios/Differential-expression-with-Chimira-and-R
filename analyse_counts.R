library(DESeq2)
library(reshape)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(plyr)


# @author: Dimitrios Vitsios & Anton Enright

# **** Input parameters ****
args=commandArgs(trailingOnly = T)


file = args[1]
CONTROL_SAMPLE = args[2]
TREATMENT_SAMPLE = args[3]

# e.g.
#file="counts_clean.txt"
#CONTROL_SAMPLE = "CTL"
#TREATMENT_SAMPLE = "KO"

# In this case, the column names in counts_clean.txt should be:
#   CTL1	CTL2	KO1	KO2

# **************************




hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
window=function(x){
  if(is.null(pdf)){
    quartz()
  } else {
    par(mfrow=c(1,1))
  }
}
let7_family_str = "let-7|mir-98-"




# *** main analysis ***
allcounts=read.table(file,header=TRUE)
rownames(allcounts) = allcounts$files
allcounts$files = NULL
allcounts=allcounts[-nrow(allcounts),]



conds=gsub("\\d+","", colnames(allcounts))


samplelevels=as.factor(conds)
samples = levels(as.factor(conds))
no_sample_levels = length(samples)
no_samples =length(conds)

samplecolors=rainbow(no_sample_levels)
names(samplecolors) = samples


coldata = as.data.frame(t(t(conds)))
rownames(coldata)=colnames(allcounts)
colnames(coldata)='condition'


#Make the count object, normalise, dispersion and testing.
dds <- DESeqDataSetFromMatrix(countData = allcounts, colData = coldata, design = ~ condition)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)

cds=dds
write.table(counts(cds,normalized=TRUE),file=paste(file,".normalizedcounts.txt",sep=""),sep="\t")



# basic QC plots
pdf('basic_QC.pdf')

heatmap.2(cor(log2(counts(dds,normalized=TRUE)+1)),trace="none",col=hmcol,main="Sample to Sample\nCorrelation", margins=c(11,11), cexRow=1, cexCol=1, cex=0.4) 
window()
pca <- princomp(counts(dds,normalized=T))
plot(pca$loadings, col=as.factor(conds),  pch=19, cex=2, main="Sample to Sample PCA")
text(pca$loadings, as.vector(colnames(allcounts)), pos=3, cex=0.8)


# Plot the total number of reads per sample before and after normalisation
plot(1, type="n", axes=F, xlab="", ylab="")
legend("center",conds,fill=as.factor(conds),cex=0.6,horiz=TRUE)
barplot(colSums(counts(cds, normalized=F)), col=as.factor(conds), las=2,cex.names=0.4,main="Pre Normalised Counts")
barplot(colSums(counts(cds, normalized=T)), col=as.factor(conds), las=2,cex.names=0.4,main="Post Normalised Counts")

# Calculate the mean normalised counts for each Tissue
meanCounts <- t(ddply(data.frame(t(counts(cds, normalized=T))), ~ conds, colMeans))
colnames(meanCounts) <- meanCounts["conds",]
meanCounts <- meanCounts[rownames(meanCounts) != "conds" ,]
meanCounts <- data.frame(row.names=rownames(meanCounts), apply(meanCounts,2,as.numeric))


# Estimate the gene-wise dispersion and plot it against the mean expression
par(mfrow=c(1,1))
plotDispEsts(cds)
dev.off()





# ************** ************** ************** ************** ************** *******
# **************** get significantly differentially expressed miRNAs ***************
# ************** ************** ************** ************** ************** *******
get_samples_diff_expr <- function(sample1, sample2, pval_thres){
  
  diff_allcounts = allcounts[ , grepl(sample1, colnames(allcounts)) | grepl(sample2, colnames(allcounts))]  
  diff_conds=gsub("\\d","",colnames(diff_allcounts))
  
  
  
  diff_samplelevels=as.factor(diff_conds)
  diff_samples = levels(as.factor(diff_conds))
  diff_no_sample_levels = length(diff_samples)
  diff_no_samples =length(diff_conds)
  
  
  diff_samplecolors=rainbow(diff_no_sample_levels)
  names(diff_samplecolors) = diff_samples
  
  
  
  diff_coldata = as.data.frame(t(t(diff_conds)))
  rownames(diff_coldata)=colnames(diff_allcounts)
  colnames(diff_coldata)='condition'
  
  
  #Make the count object, normalise, dispersion and testing.
  diff_dds <- DESeqDataSetFromMatrix(countData = diff_allcounts, colData = diff_coldata, design = ~ condition)
  diff_dds <- estimateSizeFactors(diff_dds)
  diff_dds <- estimateDispersions(diff_dds)
  diff_dds <- nbinomWaldTest(diff_dds)
  
  
  diff_dds = DESeq(diff_dds)
  res <- results(diff_dds)
  res <- res[order(res$padj),]
  
  
  # select the adjusted p-value for statistical signifance evaluation of differential expression.
  res_df = subset(res, select=c(log2FoldChange, padj) )
  res_df = as.data.frame(res_df)
  
  res_df[is.na(res_df)] = 1
  print(head(res_df))
  str(res_df)
  
  res_df$aux = 1
  res_signif = as.data.frame(res_df[ res_df$padj < pval_thres, ])
  res_signif$aux = NULL
  res_signif = res_signif[ abs(res_signif$log2FoldChange) >= 1,  ]
  
  
  #   print(mcols(res, use.names=TRUE))
  print("Num of significantly diff. expressed miRNAs:")
  print(nrow(res_signif))
  
  
  res_to_report_df = as.data.frame(res_df[ res_df$padj < pval_thres, ])
  res_to_report_df = res_to_report_df[ with(res_to_report_df, order(-log2FoldChange)), ]
  res_to_report_df = res_to_report_df[ abs(res_to_report_df$log2FoldChange ) >= 1, ]
  
  res_to_report_df$aux = NULL
  res_to_report_df$miRNA = rownames(res_to_report_df)
  res_to_report_df = res_to_report_df[ , c(3, 1, 2)]
  
  write.table(res_to_report_df, file=paste(sample1, "-", sample2, "_diff_expr_miRNAs.txt", sep=""), row.names=F, quote=F, sep="\t")
  
  return(res_signif)
}
# *************************************************************************
# *************************************************************************




# ___ 1 ___
# **** Scatterplot for miRNAs overall expression across two samples ****
norm_df = counts(cds,normalized=TRUE)

sample_general_ids = c('sample1', 'sample2')

norm_df=transform(norm_df, sample1 = eval(parse(text=paste(CONTROL_SAMPLE, 1, sep='')))  + eval(parse(text=paste(CONTROL_SAMPLE, 2, sep='')))) 
norm_df=transform(norm_df, sample2 = eval(parse(text=paste(TREATMENT_SAMPLE, 1, sep=''))) + eval(parse(text=paste(TREATMENT_SAMPLE, 2, sep='')))) 

norm_df=norm_df[ , colnames(norm_df) %in% sample_general_ids]



log_norm_df = log2(norm_df + 1)
names(log_norm_df) = c(CONTROL_SAMPLE, TREATMENT_SAMPLE)
  
log_norm_df$label = rownames(log_norm_df)
log_norm_df$colour_cond = 0
log_norm_df[grepl(let7_family_str, rownames(log_norm_df)), "colour_cond"] = 100



# *********************** ----------------------- *************************
# get list of the statistically significant differentially expressed miRNAs.
ctl_treatm_signif_df = get_samples_diff_expr(CONTROL_SAMPLE, TREATMENT_SAMPLE, 0.05)

log_norm_df = merge(log_norm_df, ctl_treatm_signif_df, by=0, all=TRUE)
log_norm_df[ is.na(log_norm_df)] = -1
rownames(log_norm_df) = log_norm_df[ , "Row.names"]
log_norm_df[ , "Row.names"] = NULL


# *********************** ----------------------- *************************


pdf('Differential_Expression_analysis.pdf')


p <- ggplot(log_norm_df, aes(eval(parse(text=CONTROL_SAMPLE)), eval(parse(text=TREATMENT_SAMPLE)), colour = padj))+
        geom_abline(intercept=0, slope=1) +
        geom_smooth(method = "lm", se = TRUE, color="red") + 
        geom_point()+xlab(CONTROL_SAMPLE)+ylab(TREATMENT_SAMPLE) +
        ggtitle(paste("miRNAs overall expression in ", CONTROL_SAMPLE, " and ", TREATMENT_SAMPLE, " samples\n(log2 of normalised counts)", sep='') ) +
        theme(plot.title = element_text(hjust = 0.5)) +
        scale_x_continuous(limits=c(0, 20), expand = c(0, 0.1)) +
        scale_y_continuous(limits=c(0, 20), expand = c(0, 0.1)) + geom_text(data=log_norm_df[ grepl(let7_family_str, rownames(log_norm_df)) & log_norm_df$padj>0, ], aes(eval(parse(text=CONTROL_SAMPLE)), eval(parse(text=TREATMENT_SAMPLE)), label=label), size = 3.2, hjust=1.1)
print(p)
ggsave(plot=p,height=10,width=10,dpi=300, filename=paste("./", CONTROL_SAMPLE, "-vs-", TREATMENT_SAMPLE, "_overall_expression_change.pdf", sep=''), useDingbats=FALSE)




  

# **** let7 family members expression ****
let7_norm_df = norm_df[grepl(let7_family_str, rownames(norm_df)), ]

# scatterplot
log_let7_norm_df = log2(let7_norm_df + 1)
p <- ggplot(log_let7_norm_df, aes(log_let7_norm_df$sample1,log_let7_norm_df$sample2, colour = sample1))+
  scale_x_continuous(limits=c(0, 20), expand = c(0, 0.1)) +
  scale_y_continuous(limits=c(0, 20), expand = c(0, 0.1)) +
  geom_abline(intercept=0, slope=1)+geom_smooth(method = "lm", se = TRUE, color="red")+geom_point()+xlab(CONTROL_SAMPLE)+ylab(TREATMENT_SAMPLE) +
  ggtitle(paste("let7 miRNAs overall expression in", CONTROL_SAMPLE, "and",  TREATMENT_SAMPLE, "samples\n(log2 of normalised counts)") ) + 
  theme(plot.title = element_text(hjust = 0.5))
print(p)



# grouped barplot with ggplot
par(mai=c(2,1,1,1))
log_let7_norm_df = log_let7_norm_df[with(log_let7_norm_df, order(-sample1)),  ]

barplot(t(as.matrix(log_let7_norm_df)), beside=TRUE, las = 2, col = c("#303D98", "#F68A33"), legend = c(CONTROL_SAMPLE, TREATMENT_SAMPLE), main=paste("let7 miRNAs overall expression in", CONTROL_SAMPLE, "and",  TREATMENT_SAMPLE, "samples"), ylab="log2 of normalised counts")




# ___ 2 ___
# ***** let7 expression across all samples *****
norm_df = counts(cds,normalized=TRUE)


let7_norm_df = as.data.frame(norm_df[grepl(let7_family_str, rownames(norm_df)), ])
let7_norm_df = let7_norm_df[with(let7_norm_df, order(-eval(parse(text=paste(CONTROL_SAMPLE, 1, sep=''))), -eval(parse(text=paste(CONTROL_SAMPLE, 2, sep=''))))), ]

write.table(let7_norm_df,file="./let7_counts.txt", quote=F, sep="\t")



# scatterplot
log_let7_norm_df = log2(let7_norm_df + 1)
write.table(log_let7_norm_df,file="./let7_log2_counts.txt", quote=F, sep="\t")


colnames(log_let7_norm_df) = c(paste(CONTROL_SAMPLE,"- repl. 1"), paste(CONTROL_SAMPLE,"- repl. 2"), paste(TREATMENT_SAMPLE, "- repl. 1"), paste(TREATMENT_SAMPLE, "- repl. 2"))
# grouped barplot with ggplot
par(mai=c(2,1,1,1))
barplot(t(as.matrix(log_let7_norm_df)), beside=TRUE, las = 2, col = c("#303D98", "#303D98", "#F68A33", "#F68A33"), legend = colnames(log_let7_norm_df), main="let7 miRNAs overall expression across all sample replicates", ylab="log2 of normalised counts")




dev.off()
