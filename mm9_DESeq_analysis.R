setwd("/shell/RNASeq/github/Elsasser_2015_Nature_RNA-Seq/")
library(corrplot)

# MAPPING: with botwie2 to mm9 genome using --fast-local)
# COUNTS: with bedtools multibam using bed filed from RefSeq and RepMasker tracks downloaded from USCS)

# ------- H3.3 enriched elements as reported in Fig 1d and Extended Data Fig 1, Elsässer et. al. 2015)
H33_enriched_EFig1 <- c("REP_ETnERV-int","REP_ETnERV2-int","REP_ETnERV3-int","REP_IAP-d-int","REP_IAPEy-int","REP_IAPEY3-int","REP_IAPEz-int","REP_IAPLTR2_Mm","REP_IAPLTR3-int","REP_IAPLTR4_I","REP_LTRIS2","REP_MMERVK10C-int","REP_MMETn-int","REP_RLTR10-int","REP_RLTR1B","REP_RLTR1B-int","REP_RLTR4_MM-int","REP_RLTR6-int")

# ERV families as in RepMasker (clustered and named according to the internal element)
ERVs <- c("REP_RLTR1B-int","REP_MMERVK10C-int","REP_MuRRS4-int","REP_RLTR4_MM-int","REP_RLTR6-int","REP_MMERGLN-int","REP_MuLV-int","REP_IAPEy-int","REP_IAPEz-int","REP_IAPLTR3-int","REP_MERVL-int","REP_MERVL_2A-int","REP_ERVL-B4-int","REP_MuRRS-int","REP_RLTR10-int","REP_MURVY-int","REP_RNERVK23-int","REP_RLTR44-int","REP_MYSERV-int","REP_IAP-d-int","REP_RMER16-int","REP_MMVL30-int","REP_ETnERV2-int","REP_RMER6-int","REP_ETnERV3-int","REP_MurERV4-int","REP_ERVL-int","REP_MTB-int","REP_RLTR45-int","REP_RMER1A","REP_IAPEY3-int","REP_ETnERV-int","REP_RodERV21-int","REP_MMTV-int","REP_RLTR22_Mus","REP_MYSERV16_I-int","REP_RLTR14-int","REP_SRV_MM-int","REP_MYSERV6-int","REP_RLTR19-int","REP_MMETn-int","REP_RMER3D-int","REP_RMER17C-int","REP_MurERV4_19-int","REP_MERVK26-int","REP_RLTR42-int","REP_RMER15-int","REP_RLTR18-int","REP_RLTR22_Mur")

# ------- Read dataset descriptions and read count files
header <- read.table("header.txt", sep="\t", header = TRUE)

RepBase_counts <- read.table("RepBase_counts_2kb.txt", sep="\t", col.names = c("Chr","Start","Stop","Gene","Score","Dir",colnames(header)), header = FALSE)[,c(4,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23)]
RepBase_counts$Copies <- 1;
RepBase_aggr <- aggregate(.~Gene, data=RepBase_counts, FUN=sum)
rownames(RepBase_aggr) <- paste("REP_",RepBase_aggr$Gene, sep="");
remove("RepBase_counts")

write.table(RepBase_aggr, pipe("pbcopy"), sep="\t")

RefSeq_counts <- read.table("RefSeq_counts.txt", sep="\t",col.names = c("Chr","Start","Stop","Gene","Score","Dir","","","","","","",colnames(header)), header = FALSE)[,c(4,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29)]
RefSeq_counts$Copies <- 1;
RefSeq_aggr <- aggregate(.~Gene, data=RefSeq_counts, FUN=sum)
rownames(RefSeq_aggr) <- paste("REF_",RefSeq_aggr$Gene, sep="");
remove("RefSeq_counts")

# ------- H3.3 datasets
CountTable <- rbind(RepBase_aggr,RefSeq_aggr)
colSums(CountTable[,2:7])

# ------- filter low expressed transcripts)
CountTable <- CountTable[apply(CountTable[,2:7],1,min)>10,]
colSums(CountTable[,2:7])

pdf("H33_Correlation_Datasets.pdf")

mcor <- cor(data.matrix(CountTable[,2:7]), method="pearson")
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(mcor, method="color", col=col(100),  
         type="upper", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         diag=TRUE,
)

mcor <- cor(data.matrix(CountTable[,2:7]), method="spearman")
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(mcor, method="color", col=col(100),  
         type="upper", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         diag=TRUE,
)
dev.off()


# ------- save copy number info for later
CountTable_Copies <- CountTable$Copies
CountTable_Gene <- CountTable$Gene
CountTable$Copies <- NULL
CountTable$Gene <- NULL

# ------- define series and conditions: H3.3 KO, H3.3 KO/KD ('def') and WT/control ('ctrl') cell lines
Des <- data.frame(row.names = colnames(CountTable[,1:6]),condition = c("def","def","ctrl","def","def","ctrl"), series= c("KD","KD","KD","KO","KO","KO"))


# --------- DESeq analysis is done in two alternative ways:
# --------- 1) evaluating series independently
# --------- 2) combining all replicates

library("DESeq")
cds = newCountDataSet(CountTable[,1:6],Des)
cds = estimateSizeFactors(cds)
sizeFactors(cds)
cds = estimateDispersions(cds)
plotDispEsts(cds, main="DESeq: Per-gene dispersion estimates")
print(plotPCA(varianceStabilizingTransformation(cds), intgroup=c("condition", "series")))


# ------- 1) test for differences in series and conditions
fit1 = fitNbinomGLMs( cds, count ~ series + condition )
fit0 = fitNbinomGLMs( cds, count ~ series  )
pvalsGLM = nbinomGLMTest( fit1, fit0 )


# ------- 2) test for differences combining all replicates
out <- nbinomTest(cds, "ctrl", "def" )

# ------- summarizing results in table
summary <- data.frame(Gene=out$id)
summary$baseMean <- out$baseMean
summary$log2FoldChange <- out$log2FoldChange
summary$pval <- out$pval
summary$pvalsGLM <- pvalsGLM

summary$Copies <- CountTable_Copies
summary$REP <- grepl("REP_", summary$Gene)
rownames(out) <- out$id

write.table(summary[summary$REP,], pipe("pbcopy"), sep="\t")
write.table(summary[summary$REP,], "DESeq_rep2kb_combined.txt", sep="\t")

out$REP <- grepl("REP_", out$id)
out$Gene <- CountTable_Gene
out$Copies <- CountTable_Copies

# ------- plotting
pdf("Fig1_H33_combined_analysis.pdf")

with(out, plot(log2FoldChange, -log10(pval), pch=20, main="H3.3 combined: Volcano plot combined", xlim=c(-3,3),ylim=c(0,20),cex=0.3,col="darkgray"))
with(subset(out, REP==TRUE ), points(log2FoldChange, -log10(pval), pch=20, cex=1, col="black"))
with(subset(out, (REP==TRUE & pval<0.05)), points(log2FoldChange, -log10(pval), pch=20, cex=1, col="red"))
with(subset(out, (REP==TRUE & pval<0.05)), text(log2FoldChange, -log10(pval), labels = Gene, pos = 4, cex=0.7, col="red"))

with(out, plot(baseMean, log2FoldChange, pch=20, main="H3.3 combined: MA plot", xlim=c(10,1000000), ylim=c(-3,3), log="x", cex=0.1,col="darkgray"))
with(subset(out, REP==TRUE ), points(baseMean, log2FoldChange, pch=20, cex=1, col="black"))
with(subset(out, (REP==TRUE & pval<0.05)), points(baseMean, log2FoldChange, pch=20, cex=1, col="red"))
with(subset(out, (REP==TRUE & pval<0.05)), text(baseMean, log2FoldChange, labels = Gene, pos = 4, cex=0.7, col="red"))
abline(h=0, v=0, lty="11")

with(out, plot(baseMean, log2FoldChange, pch=20, main="H3.3 combined: MA plot", xlim=c(0,40000), ylim=c(-3,3), cex=0.1,col="darkgray"))
with(subset(out, REP==TRUE ), points(baseMean, log2FoldChange, pch=20, cex=1, col="black"))
with(subset(out, (REP==TRUE & pval<0.05)), points(baseMean, log2FoldChange, pch=20, cex=1, col="red"))
with(subset(out, (REP==TRUE & pval<0.05)), text(baseMean, log2FoldChange, labels = Gene, pos = 4, cex=0.7, col="red"))
abline(h=0, v=0, lty="11")

with(out, plot(Copies, baseMean, pch=20, main="Expression vs Copy Number (ERVs: blue, significant: red)", log="xy", ylim=c(10,1000000), cex=0.1,col="darkgray"))
with(subset(out, REP==TRUE ), points(Copies, baseMean, pch=20, cex=1, col="black"))
with(out[ERVs,], points(Copies, baseMean, pch=20, cex=1, col="blue"))
with(subset(out, (REP==TRUE & pval<0.05)), points(Copies, baseMean, pch=20, cex=1, col="red"))
with(subset(out, (REP==TRUE & pval<0.05)), text(Copies, baseMean, labels = Gene, pos = 4, cex=0.7, col="blue"))

out$perCopyBaseMean <- out$baseMean/out$Copies
with(out, plot(Copies, perCopyBaseMean, pch=20, main="Expression vs Copy Number (ERVs: blue, significant: red)", log="xy", ylim=c(0.5,10000), cex=0.1,col="darkgray"))
with(subset(out, REP==TRUE ), points(Copies, perCopyBaseMean, pch=20, cex=1, col="black"))
with(out[ERVs,], points(Copies, perCopyBaseMean, pch=20, cex=1, col="blue"))
with(subset(out, (REP==TRUE & pval<0.05)), points(Copies, perCopyBaseMean, pch=20, cex=1, col="red"))
with(subset(out, (REP==TRUE & pval<0.05)), text(Copies, perCopyBaseMean, labels = Gene, pos = 4, cex=0.7, col="red"))


with(out, plot(baseMean, log2FoldChange, pch=20, main="H3.3 KO/KD: ERVs with (blue) and without H3.3 (black)", xlim=c(10,1000000), ylim=c(-3,3), log="x", cex=0.1,col="darkgray"))
with(out[ERVs,], points(baseMean, log2FoldChange, pch=20, cex=1, col="black"))
with(out[H33_enriched_EFig1,], points(baseMean, log2FoldChange, pch=20, cex=1, col="blue"))
with(out[H33_enriched_EFig1,], text(baseMean, log2FoldChange, labels = Gene, pos = 4, cex=0.7, col="blue"))
abline(h=0, v=0, lty="11")

dev.off();dev.off()

# ------- DAXX datasets
# --------- analysis

CountTable <- rbind(RepBase_aggr,RefSeq_aggr)

# ------- filter low expressed transcripts (same as for H3.3 for consistency)
CountTable <- CountTable[apply(CountTable[,2:7],1,min)>10,]

head(CountTable[,c(8,9,13,14,17,18)])

pdf("DAXX_Correlation_Datasets.pdf")

mcor <- cor(data.matrix(CountTable[,c(8,9,13,14,17,18)]), method="pearson")
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(mcor, method="color", col=col(100),  
         type="upper", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         diag=TRUE,
)

mcor <- cor(data.matrix(CountTable[,c(8,9,13,14,17,18)]), method="spearman")
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(mcor, method="color", col=col(100),  
         type="upper", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         diag=TRUE,
)
dev.off()

# --------- DESEq
Des <- data.frame(row.names = colnames(CountTable[,c(8,9,13,14,17,18)]),condition = c("def","ctrl","def","ctrl","def","ctrl"), series= c("S1","S1","S2","S2","S3","S3"))
cds = newCountDataSet(CountTable[,c(8,9,13,14,17,18)],Des)
cds = estimateSizeFactors(cds)
sizeFactors(cds)
cds = estimateDispersions(cds, method='blind',sharingMode="fit-only",fitType="local")
plotDispEsts(cds, main="DESeq: Per-gene dispersion estimates")
print(plotPCA(varianceStabilizingTransformation(cds), intgroup=c("condition", "series")))
fit1 = fitNbinomGLMs( cds, count ~ series + condition )
fit0 = fitNbinomGLMs( cds, count ~ series  )
pvalsGLM = nbinomGLMTest( fit1, fit0 )


cds = newCountDataSet(CountTable[,c(8,9,13,14,17,18)],c("def","ctrl","def","ctrl","def","ctrl"))
cds = estimateSizeFactors(cds)
sizeFactors(cds)
cds = estimateDispersions(cds)

out <- nbinomTest(cds, "ctrl", "def" )

out$Gene <- CountTable_Gene
out$Copies <- CountTable_Copies
out$REP <- grepl("REP_", out$id)
out$pvalsGLM <- pvalsGLM
write.table(cbind(out,CountTable), pipe("pbcopy"), sep="\t", row.names=T)
summary$DAXX_log2FoldChange=out$log2FoldChange
summary$DAXX_pval=out$pvalsGLM
row.names(out) <- out$id

# ------- plotting
pdf("Fig2_DAXX_analysis_3_datasets.pdf")

with(out, plot(log2FoldChange, -log10(pval), pch=20, main="DAXX KO: Volcano plot", xlim=c(-3,3),ylim=c(0,5),cex=0.3,col="darkgray"))
with(subset(out, REP==TRUE ), points(log2FoldChange, -log10(pval), pch=20, cex=1, col="black"))
with(subset(out, (REP==TRUE & pvalsGLM<0.05)), points(log2FoldChange, -log10(pval), pch=20, cex=1, col="red"))
with(subset(out, (REP==TRUE & pvalsGLM<0.05)), text(log2FoldChange, -log10(pval), labels = Gene, pos = 4, cex=0.5, col="red"))

with(out, plot(baseMean, log2FoldChange, pch=20, main="DAXX KO: MA plot", xlim=c(10,1000000), ylim=c(-3,3), log="x", cex=0.1,col="darkgray"))
with(subset(out, REP==TRUE ), points(baseMean, log2FoldChange, pch=20, cex=1, col="black"))
with(subset(out, (REP==TRUE & pvalsGLM<0.05)), points(baseMean, log2FoldChange, pch=20, cex=1, col="red"))
#with(subset(out, (pval<0.05)), points(baseMean, log2FoldChange, pch=20, cex=0.1, col="coral"))
with(subset(out, (REP==TRUE & pvalsGLM<0.05)), text(baseMean, log2FoldChange, labels = Gene, pos = 4, cex=0.5, col="red"))
abline(h=0, v=0, lty="11")

with(out, plot(baseMean, log2FoldChange, pch=20, main="DAXXKO: ERVs with (blue) and without H3.3 (black)", xlim=c(10,1000000), ylim=c(-3,3), log="x", cex=0.1,col="darkgray"))
with(out[ERVs,], points(baseMean, log2FoldChange, pch=20, cex=1, col="black"))
with(out[H33_enriched_EFig1,], points(baseMean, log2FoldChange, pch=20, cex=1, col="blue"))
with(out[H33_enriched_EFig1,], text(baseMean, log2FoldChange, labels = Gene, pos = 4, cex=0.5, col="blue"))
abline(h=0, v=0, lty="11")


with(summary, plot(log2FoldChange, DAXX_log2FoldChange, pch=20, main="H3.3 KO vs DAXX KO: p<0.05 in H3.3 KO",cex=0.3,col="lightgray", xlim=c(-3,3), ylim=c(-3,3)))
with(subset(summary, REP==TRUE ), points(log2FoldChange, DAXX_log2FoldChange, pch=20, cex=1, col="black"))
with(subset(summary, (REP==TRUE & pval<0.05)), points(log2FoldChange, DAXX_log2FoldChange, pch=20, cex=1, col="red"))
with(subset(summary, (REP==TRUE & pval<0.05)), text(log2FoldChange, DAXX_log2FoldChange, labels = Gene, pos = 4, cex=0.7, col="red"))
abline(h=0, v=0, lty="11")

with(summary, plot(log2FoldChange, DAXX_log2FoldChange, pch=20, main="H3.3 KO vs DAXX KO: p<0.05 in DAXX KO",cex=0.3,col="lightgray", xlim=c(-3,3), ylim=c(-3,3)))
with(subset(summary, REP==TRUE ), points(log2FoldChange, DAXX_log2FoldChange, pch=20, cex=1, col="black"))
with(subset(summary, (REP==TRUE & DAXX_pval<0.05)), points(log2FoldChange, DAXX_log2FoldChange, pch=20, cex=1, col="red"))
with(subset(summary, (REP==TRUE & DAXX_pval<0.05)), text(log2FoldChange, DAXX_log2FoldChange, labels = Gene, pos = 4, cex=0.7, col="red"))
abline(h=0, v=0, lty="11")

with(summary, plot(log2FoldChange, DAXX_log2FoldChange, pch=20, main="H3.3 KO vs DAXX KO: ERVs with (blue) and without H3.3 (black)",cex=0.3,col="lightgray", xlim=c(-3,3), ylim=c(-3,3)))
with(summary[ERVs,], points(log2FoldChange, DAXX_log2FoldChange, pch=20, cex=1, col="black"))
with(summary[H33_enriched_EFig1,], points(log2FoldChange, DAXX_log2FoldChange, pch=20, cex=1, col="blue"))
with(summary[H33_enriched_EFig1,], text(log2FoldChange, DAXX_log2FoldChange, labels = Gene, pos = 4, cex=0.7, col="blue"))
abline(h=0, v=0, lty="11")

dev.off();dev.off()

# ------- ATRX datasets
head(CountTable[,c(12,14)])

# --------- analysis
cds = newCountDataSet(CountTable[,c(12,14)],c("def","ctrl"))
cds = estimateSizeFactors(cds)
sizeFactors(cds)
cds = estimateDispersions(cds, method="blind", sharingMode="fit-only")

# ------- testing without taking into account the different cell lines
out <- nbinomTest(cds, "ctrl", "def" )

out$REP <- grepl("REP_", out$id)
summary$ATRX_log2FoldChange=out$log2FoldChange

# ------- plotting
pdf("ATRX_analysis.pdf")

with(out, plot(baseMean, log2FoldChange, pch=20, main="ATRXKO MA plot combined", xlim=c(10,1000000), ylim=c(-3,3), log="x", cex=0.1,col="darkgray"))
with(subset(out, REP==TRUE ), points(baseMean, log2FoldChange, pch=20, cex=1, col="black"))
with(subset(out, (REP==TRUE & pval<0.05)), points(baseMean, log2FoldChange, pch=20, cex=1, col="red"))
#with(subset(out, (pval<0.05)), points(baseMean, log2FoldChange, pch=20, cex=0.1, col="coral"))
with(subset(out, (REP==TRUE & pval<0.05)), text(baseMean, log2FoldChange, labels = id, pos = 4, cex=0.5, col="red"))
abline(h=0, v=0, lty="11")

with(summary, plot(log2FoldChange, ATRX_log2FoldChange, pch=20, main="H3.3 KO vs ATRX KO: p<0.05 in H3.3 KO",cex=0.3,col="lightgray", xlim=c(-3,3), ylim=c(-3,3)))
with(subset(summary, REP==TRUE ), points(log2FoldChange, ATRX_log2FoldChange, pch=20, cex=1, col="black"))
with(subset(summary, (REP==TRUE & pval<0.05)), points(log2FoldChange, ATRX_log2FoldChange, pch=20, cex=1, col="red"))
with(subset(summary, (REP==TRUE & pval<0.05)), text(log2FoldChange, ATRX_log2FoldChange, labels = Gene, pos = 4, cex=0.7, col="red"))
abline(h=0, v=0, lty="11")


dev.off();dev.off()
