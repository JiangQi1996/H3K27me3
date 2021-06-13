
colfunc<-colorRampPalette(c(mycolor5[4],mycolor5[8]))
library(RColorBrewer)
library(DESeq2)
mycolor<-brewer.pal(9,"Set1")
mycolor5<-brewer.pal(9,"Blues")
#peak signal correlation in NE2/KYSE450----
broad_gene <- read.table("./Figure4A.broad.signal.txt",sep="\t",header = T)
broad_gene$diff <- broad_gene$KYSE450-broad_gene$NE2
x <- densCols(broad_gene$NE2,broad_gene$KYSE450,colramp=colorRampPalette(c("black", "white")))
broad_gene$dens <- col2rgb(x)[1,] + 1L
cols <-  colorRampPalette(rev(brewer.pal(11,'RdYlBu')))(256)
broad_gene$col <- cols[broad_gene$dens]
plot(
  KYSE450 ~ NE2,
  data = broad_gene[order(broad_gene$dens), ],
  ylim = c(0,30),
  xlim = c(0,30),
  xlab = expression(bold('H3K27me3 signal (NE2)')),
  ylab = expression(bold('H3K27me3 signal (KYSE450)')),
  pch = 20,
  cex.lab = 1.5,
  font.lab = 4,
  col=col,
  axes = F
)
abline(a=1,b=1,lty=2,col="black",lwd=2)
abline(a=-1,b=1,lty=2,col="black",lwd=2)
axis(
  2,
  lwd = 3,
  lwd.ticks = 3,
  cex.axis = 1.5,
  font.axis = 1
)
axis(
  1,
  lwd = 3,
  lwd.ticks = 3,
  cex.axis = 1.5,
  font.axis = 1
)


narrow_gene <- read.table("./Figure4A.narrow.signal.txt",sep="\t",header = T)
narrow_gene$diff <- narrow_gene$KYSE450-narrow_gene$NE2
x <- densCols(narrow_gene$NE2,narrow_gene$KYSE450,colramp=colorRampPalette(c("black", "white")))
narrow_gene$dens <- col2rgb(x)[1,] + 1L
cols <-  colorRampPalette(rev(brewer.pal(11,'RdYlBu')))(256)
narrow_gene$col <- cols[narrow_gene$dens]
plot(
  KYSE450 ~ NE2,
  data = narrow_gene[order(narrow_gene$dens), ],
  ylim = c(0,15),
  xlim = c(0,15),
  xlab = expression(bold('H3K27me3 signal (NE2)')),
  ylab = expression(bold('H3K27me3 signal (KYSE510)')),
  pch = 20,
  cex.lab = 1.5,
  font.lab = 4,
  col=col,
  axes = F
)
abline(a=1,b=1,lty=2,col="black",lwd=2)
abline(a=-1,b=1,lty=2,col="black",lwd=2)
axis(
  2,
  lwd = 3,
  lwd.ticks = 3,
  cex.axis = 1.5,
  font.axis = 1
)
axis(
  1,
  lwd = 3,
  lwd.ticks = 3,
  cex.axis = 1.5,
  font.axis = 1
)

random_gene <- read.table("./Figure4A.random.signal.txt",sep="\t",header = T)
random_gene$diff <- random_gene$KYSE450-random_gene$NE2
x <- densCols(random_gene$NE2,random_gene$KYSE450,colramp=colorRampPalette(c("black", "white")))
random_gene$dens <- col2rgb(x)[1,] + 1L
cols <-  colorRampPalette(rev(brewer.pal(11,'RdYlBu')))(256)
random_gene$col <- cols[random_gene$dens]
plot(
  KYSE450 ~ NE2,
  data = random_gene[order(random_gene$dens), ],
  ylim = c(0,20),
  xlim = c(0,20),
  xlab = expression(bold('H3K27me3 signal (NE2)')),
  ylab = expression(bold('H3K27me3 signal (KYSE510)')),
  pch = 20,
  cex.lab = 1.5,
  font.lab = 4,
  col=col,
  axes = F
)
abline(a=1,b=1,lty=2,col="black",lwd=2)
abline(a=-1,b=1,lty=2,col="black",lwd=2)
axis(
  2,
  lwd = 3,
  lwd.ticks = 3,
  cex.axis = 1.5,
  font.axis = 1
)
axis(
  1,
  lwd = 3,
  lwd.ticks = 3,
  cex.axis = 1.5,
  font.axis = 1
)

#DESeq2----
ensg2gene <- read.table("./ensg2gene",stringsAsFactors = F)
gene_pcg <- read.table("I:/H3K27me3/Homo_sapiens.GRCh38.95.gtf.pcg",stringsAsFactors = F,sep="\t")
dir(path = "I:/H3K27me3/",
    pattern = "*htseq.count",
    full.names = T) -> ESCC_gene_htseq
dir(path = "I:/H3K27me3/",
    pattern = "*htseq.count",
    full.names = F) -> sample
sample <- sub(".htseq.count", "", sample)
lapply(
  ESCC_gene_htseq,
  read.table,
  header = F,
  col.names = c("gene","count"),
  quote = "",
  sep = "\t",
  as.is = T
) -> ESCC_gene_htseq_all
total_row1 <- unique(unlist(lapply(ESCC_gene_htseq_all, function(x) {
  return(x$gene)
})))
raw_df1 = data.frame(gene = total_row1)
for (file in ESCC_gene_htseq_all) {
  raw_df1 = cbind(raw_df1, file$count)
}
colnames(raw_df1) <- c("gene", sample)
ESCC_gene_count <- raw_df1[grep("ENSG", raw_df1$gene), ]
dim(ESCC_gene_count)
#[1] 58735     7
ESCC_gene_count <- ESCC_gene_count[, 2:7]
rownames(ESCC_gene_count) <- raw_df1[grep("ENSG", raw_df1$gene), 1]
dim(ESCC_gene_count)
#[1] 58735    6
##过滤count都很低的基因
ESCC_gene_count1 <- ESCC_gene_count
percentile <- function(x) {
  length(x[x >= 10]) / 6
}
ESCC_gene_count1$percentile <- apply(ESCC_gene_count1, 1, percentile)
ESCC_gene_count1 <-
  ESCC_gene_count1[which(ESCC_gene_count1$percentile >= 0.05), 1:6]

ESCC_gene_count1 <- ESCC_gene_count1[which(rownames(ESCC_gene_count1) != "ENSG00000277483"),]
dim(ESCC_gene_count1)
#[1] 23974     6
#Normalization----
colData_all <- data.frame(rep(c("tumor_450","tumor_510","normal"),each=2),"pair_end")
colnames(colData_all) <- c("condition", "type")
rownames(colData_all) <- c("R18043758LR01-KYSE450-1","R18043759LR01-KYSE450-2",
                           "R18056087LR01-KYSE510-1","R18056088LR01-KYSE510-2",
                           "R18056081LR01-NE2-1","R18056082LR01-NE2-2")
ESCC_gene_count1 <- ESCC_gene_count1[,order(match(colnames(ESCC_gene_count1),rownames(colData_all)))]
dds <- DESeqDataSetFromMatrix(countData = ESCC_gene_count1,
                         colData = colData_all,
                         design = ~ condition)
dds <- DESeq(dds)
#get transformed dataset
vsd <- vst(dds, blind=T)
rld <- rlog(dds, blind=T)
head(assay(vsd), 5)
head(assay(rld), 5)
dcounts_gene_deseq <- counts(dds, normalized = TRUE)
ESCC_gene_deseq <- as.data.frame(log2(dcounts_gene_deseq+1))
head(ESCC_gene_deseq,5)
dim(ESCC_gene_deseq)
#[1] 23974   6
#450
DGE_tumor_450 <-
  results(
    dds,
    contrast = c("condition", "tumor_450", "normal"),
    alpha = 0.05
  )
summary(DGE_tumor_450)
plotCounts(dds, gene=which.min(DGE_tumor_450$padj), intgroup="condition")
DGE_tumor_450 <-
  data.frame(
    gene = rownames(DGE_tumor_450),
    LFC = DGE_tumor_450$log2FoldChange,
    padj = DGE_tumor_450$padj
  )

DGE_tumor_450 <- merge(DGE_tumor_450,ensg2gene,by.x="gene",by.y="V1")
DGE_tumor_450_sig <- DGE_tumor_450[which(abs(DGE_tumor_450$LFC) >=1 & DGE_tumor_450$padj<=0.05),]

#510
DGE_tumor_510 <-
  results(
    dds,
    contrast = c("condition", "tumor_510", "normal"),
    alpha = 0.05,
    cooksCutoff = F
  )
DGE_tumor_510 <-
  data.frame(
    gene = rownames(DGE_tumor_510),
    LFC = DGE_tumor_510$log2FoldChange,
    padj = DGE_tumor_510$padj
  )
DGE_tumor_510 <- merge(DGE_tumor_510,ensg2gene,by.x="gene",by.y="V1")
DGE_tumor_510_sig <- DGE_tumor_510[which(abs(DGE_tumor_510$LFC) >=1 & DGE_tumor_510$padj<=0.05),]
DGE_tumor <- merge(DGE_tumor_450,DGE_tumor_510,by="gene")

##correlation between peak signal and mRNA expression----
DGE_tumor_broad1 <- merge(DGE_tumor,broad_gene,by.x="gene",by.y="geneid")
DGE_tumor_broad1 <- unique(DGE_tumor_broad1[,c("gene","LFC.x","padj.x","symbol","NE2","KYSE450","diff")])
a <- cor.test(DGE_tumor_broad1$LFC.x,DGE_tumor_broad1$diff,method = "spearman")
DGE_tumor_broad2 <- DGE_tumor_broad1[which(DGE_tumor_broad1$diff <= -1 & DGE_tumor_broad1$LFC.x >= 1 & DGE_tumor_broad1$padj.x <= 0.05),]
DGE_tumor_broad3 <- DGE_tumor_broad1[which(DGE_tumor_broad1$diff>=1 & DGE_tumor_broad1$LFC.x <= -1 & DGE_tumor_broad1$padj.x <= 0.05),]
plot(DGE_tumor_broad1$diff,DGE_tumor_broad1$LFC.x)
plot(
  DGE_tumor_broad1$diff,DGE_tumor_broad1$LFC.x,
  ylim = c(-15,15),
  xlim = c(-25,15),
  xlab = expression(bold('H3K27me3 signal change')),
  ylab = expression(bold('mRNA fold change')),
  pch = 20,
  col ="grey",
  cex.lab = 1.5,
  font.lab = 4,
  axes = F
)
points(DGE_tumor_broad2$diff,DGE_tumor_broad2$LFC.x,col=mycolor[1],  pch = 20)
points(DGE_tumor_broad3$diff,DGE_tumor_broad3$LFC.x,col=mycolor[2],  pch = 20)
abline(lm(DGE_tumor_broad1$LFC.x~DGE_tumor_broad1$diff),col="black")
axis(
  2,
  lwd = 3,
  lwd.ticks = 3,
  cex.axis = 1.5,
  font.axis = 1
)
axis(
  1,
  lwd = 3,
  lwd.ticks = 3,
  cex.axis = 1.5,
  font.axis = 1
)


DGE_tumor_narrow1 <- merge(DGE_tumor,narrow_gene,by.x="gene",by.y="geneid")
DGE_tumor_narrow1 <- unique(DGE_tumor_narrow1[,c("gene","LFC.x","padj.x","symbol","NE2","KYSE450","diff")])
a <- cor.test(DGE_tumor_narrow1$LFC.x,DGE_tumor_narrow1$diff,method = "spearman")
DGE_tumor_narrow2 <- DGE_tumor_narrow1[which(DGE_tumor_narrow1$diff <= -1 & DGE_tumor_narrow1$LFC.x >= 1 & DGE_tumor_narrow1$padj.x <= 0.05),]
DGE_tumor_narrow3 <- DGE_tumor_narrow1[which(DGE_tumor_narrow1$diff>=1 & DGE_tumor_narrow1$LFC.x <= -1 & DGE_tumor_narrow1$padj.x <= 0.05),]
plot(
  DGE_tumor_narrow1$diff,DGE_tumor_narrow1$LFC.x,
  ylim = c(-10,15),
  xlim = c(-10,15),
  xlab = expression(bold('H3K27me3 signal change')),
  ylab = expression(bold('mRNA fold change')),
  pch = 20,
  col ="grey",
  cex.lab = 1.5,
  font.lab = 4,
  axes = F
)
points(DGE_tumor_narrow2$diff,DGE_tumor_narrow2$LFC.x,col=mycolor[1],  pch = 20)
points(DGE_tumor_narrow3$diff,DGE_tumor_narrow3$LFC.x,col=mycolor[2],  pch = 20)
abline(lm(DGE_tumor_narrow1$LFC.x~DGE_tumor_narrow1$diff),col="black")
axis(
  2,
  lwd = 3,
  lwd.ticks = 3,
  cex.axis = 1.5,
  font.axis = 1
)
axis(
  1,
  lwd = 3,
  lwd.ticks = 3,
  cex.axis = 1.5,
  font.axis = 1
)

DGE_tumor_random1 <- merge(DGE_tumor,random_gene,by.x="gene",by.y="geneid")
DGE_tumor_random1 <- unique(DGE_tumor_random1[,c("gene","LFC.x","padj.x","symbol","NE2","KYSE450","diff")])
a <- cor.test(DGE_tumor_random1$LFC.x,DGE_tumor_random1$diff,method = "spearman")
DGE_tumor_random2 <- DGE_tumor_random1[which(DGE_tumor_random1$diff <= -1 & DGE_tumor_random1$LFC.x >= 1 & DGE_tumor_random1$padj.x <= 0.05),]
DGE_tumor_random3 <- DGE_tumor_random1[which(DGE_tumor_random1$diff>=1 & DGE_tumor_random1$LFC.x <= -1 & DGE_tumor_random1$padj.x <= 0.05),]
plot(DGE_tumor_random1$diff,DGE_tumor_random1$LFC.x)
plot(
  DGE_tumor_random1$diff,DGE_tumor_random1$LFC.x,
  ylim = c(-15,15),
  xlim = c(-10,10),
  xlab = expression(bold('H3K27me3 signal change')),
  ylab = expression(bold('mRNA fold change')),
  pch = 20,
  col ="grey",
  cex.lab = 1.5,
  font.lab = 4,
  axes = F
)
points(DGE_tumor_random2$diff,DGE_tumor_random2$LFC.x,col=mycolor[1],  pch = 20)
points(DGE_tumor_random3$diff,DGE_tumor_random3$LFC.x,col=mycolor[2],  pch = 20)
abline(lm(DGE_tumor_random1$LFC.x~DGE_tumor_random1$diff),col="black")
axis(
  2,
  lwd = 3,
  lwd.ticks = 3,
  cex.axis = 1.5,
  font.axis = 1
)
axis(
  1,
  lwd = 3,
  lwd.ticks = 3,
  cex.axis = 1.5,
  font.axis = 1
)

#TCGA----
gene_pcg <- read.table("./Homo_sapiens.GRCh38.95.gtf.pcg",stringsAsFactors = F,sep="\t")
ESCA <- read.table("./combn_106.DEG.results.txt")
ESCA$gene <- rownames(ESCA)
ESCA <- merge(ESCA,gene_pcg,by.x="gene",by.y="V5")
ESCA_broad1 <- merge(ESCA,broad_gene,by.x="gene",by.y="geneid")
ESCA_broad1 <- unique(ESCA_broad1[,c("gene","log2FoldChange","padj","symbol","NE2","KYSE510","diff")])
ESCA_broad1 <- ESCA_broad1[which(ESCA_broad1$gene %in% DGE_tumor_broad1$gene),]
cor.test(ESCA_broad1$log2FoldChange,ESCA_broad1$diff,method = "spearman")
ESCA_broad2 <- ESCA_broad1[which(ESCA_broad1$gene %in% DGE_tumor_broad2$gene),]
plot(ESCA_broad1$diff,ESCA_broad1$log2FoldChange)
plot(
  ESCA_broad1$diff,ESCA_broad1$log2FoldChange,
  ylim = c(-10,10),
  xlim = c(-25,10),
  xlab = expression(bold('H3K27me3 signal change')),
  ylab = expression(bold('mRNA fold change')),
  pch = 20,
  col ="black",
  cex.lab = 1.5,
  font.lab = 2,
  axes = F
)
points(ESCA_broad2$diff,ESCA_broad2$log2FoldChange,col=mycolor[1],  pch = 20)
abline(lm(ESCA_broad2$log2FoldChange~ESCA_broad2$diff),col=mycolor[1])
axis(
  2,
  lwd = 3,
  lwd.ticks = 3,
  cex.axis = 1.5,
  font.axis = 1
)
axis(
  1,
  lwd = 3,
  lwd.ticks = 3,
  cex.axis = 1.5,
  font.axis = 1
)

DGE_tumor_broad1 <- merge(DGE_tumor,broad_gene,by.x="gene",by.y="geneid")
DGE_tumor_broad1 <- unique(DGE_tumor_broad1[,c("gene","LFC.x","padj.x","symbol","NE2","KYSE450","diff")])
cor.test(DGE_tumor_broad1$LFC.x,DGE_tumor_broad1$diff,method = "spearman")
DGE_tumor_broad2 <- DGE_tumor_broad1[which(DGE_tumor_broad1$diff <= -1 & DGE_tumor_broad1$LFC.x >= 1 & DGE_tumor_broad1$padj.x <= 0.05),]

DGE_tumor_narrow1 <- merge(DGE_tumor,narrow_gene,by.x="gene",by.y="geneid")
DGE_tumor_narrow1 <- unique(DGE_tumor_narrow1[,c("gene","LFC.x","padj.x","symbol","NE2","KYSE450","diff")])
cor.test(DGE_tumor_narrow1$LFC.x,DGE_tumor_narrow1$diff,method = "spearman")
DGE_tumor_narrow2 <- DGE_tumor_narrow1[which(DGE_tumor_narrow1$diff <= -1 & DGE_tumor_narrow1$LFC.x >= 1 & DGE_tumor_narrow1$padj.x <= 0.05),]

DGE_tumor_random1 <- merge(DGE_tumor,random_gene,by.x="gene",by.y="geneid")
DGE_tumor_random1 <- unique(DGE_tumor_random1[,c("gene","LFC.x","padj.x","symbol","NE2","KYSE450","diff")])
cor.test(DGE_tumor_random1$LFC.x,DGE_tumor_random1$diff,method = "spearman")
DGE_tumor_random2 <- DGE_tumor_random1[which(DGE_tumor_random1$diff <= -1 & DGE_tumor_random1$LFC.x >= 1 & DGE_tumor_random1$padj.x <= 0.05),]


ESCA_loss1 <- ESCA[which(ESCA$gene %in% DGE_tumor_broad2$gene),]
ESCA_loss2 <- ESCA[which(ESCA$gene %in% DGE_tumor_narrow2$gene),]
ESCA_loss3 <- ESCA[which(ESCA$gene %in% DGE_tumor_random2$gene),]
par(mfrow = c(1, 1),mar=c(6,7,5,2),lwd=2)
cumplot<-function(x,colo="black",lwd=1,xlab="",ylab="",main="",xpd=F,cexmain=2,cexlab=1.8,fontlab=2,fontmain=2){x.o<-sort(x);plot(x.o,1:length(x.o)/length(x.o),xlim=c(-10,10),ylim=c(0.0,1.0),type="l",col=colo,cex.main=cexmain,cex.lab=cexlab,font.lab=fontlab,font.main=fontmain,lwd=lwd,main=main,xlab=xlab,ylab=ylab,axes=F)}
cumlines<-function(x,colo="black",lwd=1){x.o<-sort(x);lines(x.o,1:length(x.o)/length(x.o),col=colo,lwd=lwd)}
cumplot(ESCA$log2FoldChange,colo="grey48",lwd=3,main="")
cumlines(ESCA_loss1$log2FoldChange,colo=mycolor[1],lwd=3)
par(new=T)
cumlines(ESCA_loss2$log2FoldChange,colo=mycolor[2],lwd=3)
par(new=T)
cumlines(ESCA_loss3$log2FoldChange,colo=mycolor[3],lwd=3)
name1<-paste("All gene ","(19645)",sep="")
name2<-paste("Broad ","(208)",sep="")
name3<-paste("Narrow ","(106)",sep="")
name4<-paste("Control ","(106)",sep="")
legend(-10,1.05,cex=1,lty=1,legend=c(name1,name2,name3,name4),bty="n",col=c("grey48",mycolor[1:3]),text.col=c("grey48",mycolor[1:3]),text.font=2)
axis(2,c(0.0,0.2,0.4,0.6,0.8,1.0),lwd=2,lwd.ticks = 2,cex.axis=2,font.axis=2)
axis(1,c(-10,-5,0,5,10),lwd=2,lwd.ticks = 2,cex.axis=2,font.axis=2)

