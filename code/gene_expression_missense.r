mycolor<-brewer.pal(9,"Blues")
mycolor2<-brewer.pal(9,"Set1")
colfunc<-colorRampPalette(c(mycolor[4],mycolor[8]))

#width and expression----
coding_peak <- read.table("./data/ChIP/NE2.coding.peak.bed",header = T)
coding_peak <- coding_peak[which(coding_peak$width>=0.5),]
coding_peak20 <- coding_peak[1:813,]
coding_peak19 <- coding_peak[(813*1+1):(813*2),]
coding_peak18 <- coding_peak[(813*2+1):(813*3),]
coding_peak17 <- coding_peak[(813*3+1):(813*4),]
coding_peak16 <- coding_peak[(813*4+1):(813*5),]
coding_peak15 <- coding_peak[(813*5+1):(813*6),]
coding_peak14 <- coding_peak[(813*6+1):(813*7),]
coding_peak13 <- coding_peak[(813*7+1):(813*8),]
coding_peak12 <- coding_peak[(813*8+1):(813*9),]
coding_peak11 <- coding_peak[(813*9+1):(813*10),]
coding_peak10 <- coding_peak[(813*10+1):(813*11),]
coding_peak9 <- coding_peak[(813*11+1):(813*12),]
coding_peak8 <- coding_peak[(813*12+1):(813*13),]
coding_peak7 <- coding_peak[(813*13+1):(813*14),]
coding_peak6 <- coding_peak[(813*14+1):(813*15),]
coding_peak5 <- coding_peak[(813*15+1):(813*16),]
coding_peak4 <- coding_peak[(813*16+1):(813*17),]
coding_peak3 <- coding_peak[(813*17+1):(813*18),]
coding_peak2 <- coding_peak[(813*18+1):(813*19),]
coding_peak1 <- coding_peak[15429:16237,]

coding_exp <- read.table("./data/ChIP/NE2.count.gene.result",header = F)
gene_biotype_pcg <- read.table("./data/ChIP/gene_biotype_pcg.uniq.bed",header = T)
coding_exp_pcg <- coding_exp[which(coding_exp$V1 %in% gene_biotype_pcg$gene),]
coding_exp_pcg$TPM <- log2(coding_exp_pcg$V2+1)
coding_exp_pcg$RPKM <- log2(coding_exp_pcg$V3+1)
coding_exp1 <- read.table("NE2.coding.peak.anno.bed",header = T)
coding_exp1 <- merge(coding_exp1,coding_exp_pcg,by.x="id.1",by.y="V1")
coding_exp1 <- coding_exp1[which(coding_exp1$id.1 %in% gene_biotype_pcg$gene),]
coding_exp1 <- coding_exp1[which(coding_exp1$id %in% coding_peak$id),]
coding_exp1$TPM <- log2(coding_exp1$V2+1)
coding_exp1$RPKM <- log2(coding_exp1$V3+1)

#with peak
coding_exp_pcg_p1 <- unique(coding_exp1[which(coding_exp1$id %in% coding_peak1$id),c(1,16:17)])
coding_exp_pcg_p2 <- unique(coding_exp1[which(coding_exp1$id %in% coding_peak2$id),c(1,16:17)])
coding_exp_pcg_p3 <- unique(coding_exp1[which(coding_exp1$id %in% coding_peak3$id),c(1,16:17)])
coding_exp_pcg_p4 <- unique(coding_exp1[which(coding_exp1$id %in% coding_peak4$id),c(1,16:17)])
coding_exp_pcg_p5 <- unique(coding_exp1[which(coding_exp1$id %in% coding_peak5$id),c(1,16:17)])
coding_exp_pcg_p6 <- unique(coding_exp1[which(coding_exp1$id %in% coding_peak6$id),c(1,16:17)])
coding_exp_pcg_p7 <- unique(coding_exp1[which(coding_exp1$id %in% coding_peak7$id),c(1,16:17)])
coding_exp_pcg_p8 <- unique(coding_exp1[which(coding_exp1$id %in% coding_peak8$id),c(1,16:17)])
coding_exp_pcg_p9 <- unique(coding_exp1[which(coding_exp1$id %in% coding_peak9$id),c(1,16:17)])
coding_exp_pcg_p10 <- unique(coding_exp1[which(coding_exp1$id %in% coding_peak10$id),c(1,16:17)])
coding_exp_pcg_p11 <- unique(coding_exp1[which(coding_exp1$id %in% coding_peak11$id),c(1,16:17)])
coding_exp_pcg_p12 <- unique(coding_exp1[which(coding_exp1$id %in% coding_peak12$id),c(1,16:17)])
coding_exp_pcg_p13 <- unique(coding_exp1[which(coding_exp1$id %in% coding_peak13$id),c(1,16:17)])
coding_exp_pcg_p14 <- unique(coding_exp1[which(coding_exp1$id %in% coding_peak14$id),c(1,16:17)])
coding_exp_pcg_p15 <- unique(coding_exp1[which(coding_exp1$id %in% coding_peak15$id),c(1,16:17)])
coding_exp_pcg_p16 <- unique(coding_exp1[which(coding_exp1$id %in% coding_peak16$id),c(1,16:17)])
coding_exp_pcg_p17 <- unique(coding_exp1[which(coding_exp1$id %in% coding_peak17$id),c(1,16:17)])
coding_exp_pcg_p18 <- unique(coding_exp1[which(coding_exp1$id %in% coding_peak18$id),c(1,16:17)])
coding_exp_pcg_p19 <- unique(coding_exp1[which(coding_exp1$id %in% coding_peak19$id),c(1,16:17)])
coding_exp_pcg_p20 <- unique(coding_exp1[which(coding_exp1$id %in% coding_peak20$id),c(1,16:17)])
coding_exp_pcg_other <- unique(coding_exp1[which(coding_exp1$id %in% setdiff(coding_exp1$id,coding_peak20$id)),c(1,15:17)])

par(new=T)
pdf(file="mRNA with H3K27me3.pdf")
boxplot(
  coding_exp_pcg_p1$TPM,
  coding_exp_pcg_p2$TPM,
  coding_exp_pcg_p3$TPM,
  coding_exp_pcg_p4$TPM,
  coding_exp_pcg_p5$TPM,
  coding_exp_pcg_p6$TPM,
  coding_exp_pcg_p7$TPM,
  coding_exp_pcg_p8$TPM,
  coding_exp_pcg_p9$TPM,
  coding_exp_pcg_p10$TPM,
  coding_exp_pcg_p11$TPM,
  coding_exp_pcg_p12$TPM,
  coding_exp_pcg_p13$TPM,
  coding_exp_pcg_p14$TPM,
  coding_exp_pcg_p15$TPM,
  coding_exp_pcg_p16$TPM,
  coding_exp_pcg_p17$TPM,
  coding_exp_pcg_p18$TPM,
  coding_exp_pcg_p19$TPM,
  coding_exp_pcg_p20$TPM,
  ylab = expression(bold('mRNA level (log2TPM)')),
  ylim = c(0,1),
  xlim = c(0,20),
  border = c(rep(mycolor[5],19),mycolor[9]),
  cex.lab = 1.5,
  font.lab = 2,
  outline = F,
  axes = F
)
axis(
  2,
  lwd = 3,
  lwd.ticks = 3,
  cex.axis = 1.5,
  font.axis = 2
)
dev.off()
#wilcox.test(coding_exp_pcg_p1$TPM,coding_exp_pcg_p20$TPM)
wilcox.test(coding_exp_pcg_p19$TPM,coding_exp_pcg_p20$TPM) #pvalue0.0008448

#essentiality
essentiality <- read.csv("./data/ChIP/nrg.2017.75-s2.txt",sep="\t",header = T)
essentiality_p1 <- essentiality[which(essentiality$X %in% coding_exp_pcg_p1$id.1),]
essentiality_p2 <- essentiality[which(essentiality$X %in% coding_exp_pcg_p2$id.1),]
essentiality_p3 <- essentiality[which(essentiality$X %in% coding_exp_pcg_p3$id.1),]
essentiality_p4 <- essentiality[which(essentiality$X %in% coding_exp_pcg_p4$id.1),]
essentiality_p5 <- essentiality[which(essentiality$X %in% coding_exp_pcg_p5$id.1),]
essentiality_p6 <- essentiality[which(essentiality$X %in% coding_exp_pcg_p6$id.1),]
essentiality_p7 <- essentiality[which(essentiality$X %in% coding_exp_pcg_p7$id.1),]
essentiality_p8 <- essentiality[which(essentiality$X %in% coding_exp_pcg_p8$id.1),]
essentiality_p9 <- essentiality[which(essentiality$X %in% coding_exp_pcg_p9$id.1),]
essentiality_p10 <- essentiality[which(essentiality$X %in% coding_exp_pcg_p10$id.1),]
essentiality_p11 <- essentiality[which(essentiality$X %in% coding_exp_pcg_p11$id.1),]
essentiality_p12 <- essentiality[which(essentiality$X %in% coding_exp_pcg_p12$id.1),]
essentiality_p13 <- essentiality[which(essentiality$X %in% coding_exp_pcg_p13$id.1),]
essentiality_p14 <- essentiality[which(essentiality$X %in% coding_exp_pcg_p14$id.1),]
essentiality_p15 <- essentiality[which(essentiality$X %in% coding_exp_pcg_p15$id.1),]
essentiality_p16 <- essentiality[which(essentiality$X %in% coding_exp_pcg_p16$id.1),]
essentiality_p17 <- essentiality[which(essentiality$X %in% coding_exp_pcg_p17$id.1),]
essentiality_p18 <- essentiality[which(essentiality$X %in% coding_exp_pcg_p18$id.1),]
essentiality_p19 <- essentiality[which(essentiality$X %in% coding_exp_pcg_p19$id.1),]
essentiality_p20 <- essentiality[which(essentiality$X %in% coding_exp_pcg_p20$id.1),]
pdf(file="missensae.Z.pdf")
boxplot(
  essentiality_p1$missense.Z,
  essentiality_p2$missense.Z,
  essentiality_p3$missense.Z,
  essentiality_p4$missense.Z,
  essentiality_p5$missense.Z,
  essentiality_p6$missense.Z,
  essentiality_p7$missense.Z,
  essentiality_p8$missense.Z,
  essentiality_p9$missense.Z,
  essentiality_p10$missense.Z,
  essentiality_p11$missense.Z,
  essentiality_p12$missense.Z,
  essentiality_p13$missense.Z,
  essentiality_p14$missense.Z,
  essentiality_p15$missense.Z,
  essentiality_p16$missense.Z,
  essentiality_p17$missense.Z,
  essentiality_p18$missense.Z,
  essentiality_p19$missense.Z,
  essentiality_p20$missense.Z,
  ylab = expression(bold('missense Z score')),
  border = c(rep(mycolor[5],19),mycolor[9]),
  cex.lab = 1.5,
  font.lab = 2,
  outline = F,
  axes = F
)
axis(
  2,
  lwd = 3,
  lwd.ticks = 3,
  cex.axis = 1.5,
  font.axis = 2
)
abline(h = 1, lty=2, col = 'red')
dev.off()
wilcox.test(essentiality_p20$missense.Z,essentiality_p19$missense.Z) #P-value = 0.0008448

#dN/dS
mouse <- read.csv("./data/ChIP/mouse.txt",header = T,sep = ",")
mouse1 <- mouse
mouse1$w <- mouse1$dN.with.Mouse/mouse1$dS.with.Mouse
mouse1 <- mouse1[which(mouse1$Mouse.homology.type=="ortholog_one2one"),]
mouse1_p1 <- mouse1[which(mouse1$Gene.stable.ID %in% coding_exp_pcg_p1$id.1),]
mouse1_p2 <- mouse1[which(mouse1$Gene.stable.ID %in% coding_exp_pcg_p2$id.1),]
mouse1_p3 <- mouse1[which(mouse1$Gene.stable.ID %in% coding_exp_pcg_p3$id.1),]
mouse1_p4 <- mouse1[which(mouse1$Gene.stable.ID %in% coding_exp_pcg_p4$id.1),]
mouse1_p5 <- mouse1[which(mouse1$Gene.stable.ID %in% coding_exp_pcg_p5$id.1),]
mouse1_p6 <- mouse1[which(mouse1$Gene.stable.ID %in% coding_exp_pcg_p6$id.1),]
mouse1_p7 <- mouse1[which(mouse1$Gene.stable.ID %in% coding_exp_pcg_p7$id.1),]
mouse1_p8 <- mouse1[which(mouse1$Gene.stable.ID %in% coding_exp_pcg_p8$id.1),]
mouse1_p9 <- mouse1[which(mouse1$Gene.stable.ID %in% coding_exp_pcg_p9$id.1),]
mouse1_p10 <- mouse1[which(mouse1$Gene.stable.ID %in% coding_exp_pcg_p10$id.1),]
mouse1_p11 <- mouse1[which(mouse1$Gene.stable.ID %in% coding_exp_pcg_p11$id.1),]
mouse1_p12 <- mouse1[which(mouse1$Gene.stable.ID %in% coding_exp_pcg_p12$id.1),]
mouse1_p13 <- mouse1[which(mouse1$Gene.stable.ID %in% coding_exp_pcg_p13$id.1),]
mouse1_p14 <- mouse1[which(mouse1$Gene.stable.ID %in% coding_exp_pcg_p14$id.1),]
mouse1_p15 <- mouse1[which(mouse1$Gene.stable.ID %in% coding_exp_pcg_p15$id.1),]
mouse1_p16 <- mouse1[which(mouse1$Gene.stable.ID %in% coding_exp_pcg_p16$id.1),]
mouse1_p17 <- mouse1[which(mouse1$Gene.stable.ID %in% coding_exp_pcg_p17$id.1),]
mouse1_p18 <- mouse1[which(mouse1$Gene.stable.ID %in% coding_exp_pcg_p18$id.1),]
mouse1_p19 <- mouse1[which(mouse1$Gene.stable.ID %in% coding_exp_pcg_p19$id.1),]
mouse1_p20 <- mouse1[which(mouse1$Gene.stable.ID %in% coding_exp_pcg_p20$id.1),]
pdf(file="dN_dS.pdf")
boxplot(
  mouse1_p1$w,
  mouse1_p2$w,
  mouse1_p3$w,
  mouse1_p4$w,
  mouse1_p5$w,
  mouse1_p6$w,
  mouse1_p7$w,
  mouse1_p8$w,
  mouse1_p9$w,
  mouse1_p10$w,
  mouse1_p11$w,
  mouse1_p12$w,
  mouse1_p13$w,
  mouse1_p14$w,
  mouse1_p15$w,
  mouse1_p16$w,
  mouse1_p17$w,
  mouse1_p18$w,
  mouse1_p19$w,
  mouse1_p20$w,
  ylab = expression(bold('dN/dS')),
  border = c(rep(mycolor[5],19),mycolor[9]),
  cex.lab = 1.5,
  font.lab = 2,
  outline = F,
  axes = F
)
axis(
  2,
  lwd = 3,
  lwd.ticks = 3,
  cex.axis = 1.5,
  font.axis = 2
)
abline(h = 0.08, lty=2, col = 'red')
dev.off()
wilcox.test(mouse1_p20$w,mouse1_p19$w)  # P-value = 0.0397