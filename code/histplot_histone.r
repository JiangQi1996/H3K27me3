#hist plot for NE2 histone peak distribution 
#H3K27me3
library(RColorBrewer)
NE2_K27me3 <- read.table("./data/ChIP/NE2H3K27me3.length.peak.bed",sep = "\t",header = F)
NE2_K27me3_1 <- NE2_K27me3[which(NE2_K27me3$V4 > 0.5 ),]
NE2_K27me3_2 <- NE2_K27me3_1[which(NE2_K27me3_1$V4 <= 80 ),]
x <- NE2_K27me3_2$V4
mycolor<-brewer.pal(9,"Set1")
pdf(file="Distribution of H3K27me3 ChIP-seq peaks in NE2.pdf")
h<-hist(x, breaks = 400, col=mycolor[2], xlab="H3K27me3 breadth (kb)",
        main="Distribution of H3K27me3 ChIP-seq peaks in NE2",border = NA,  font.lab = 2, cex.lab = 1.2)
dev.off()

#H3K27ac
NE2_K27ac <- read.table("./data/ChIP/NE2H3K27ac.length.peak.bed",sep = "\t",header = F)
NE2_K27ac_2 <- NE2_K27ac[which(NE2_K27ac$V4 <= 20 ),]
x <- NE2_K27ac_2$V4
mycolor<-brewer.pal(9,"Set1")
pdf(file="Distribution of H3K27ac ChIP-seq peaks in NE21.pdf")
h<-hist(x, breaks = 100, col=mycolor[2], xlab="H3K27ac breadth (kb)",
        main="Distribution of H3K27ac ChIP-seq peaks in NE2",border = NA,  font.lab = 2, cex.lab = 1.2)
dev.off()

#H3K4me3
NE2_K4me3 <- read.table("./data/ChIP/NE2H3K4me3.length.peak.bed",sep = "\t",header = F)
NE2_K4me3_2 <- NE2_K4me3[which(NE2_K4me3$V4 <= 25 ),]
x <- NE2_K4me3_2$V4
mycolor<-brewer.pal(9,"Set1")
pdf(file="Distribution of H3K4me3 ChIP-seq peaks in NE2.pdf")
h<-hist(x, breaks = 100, col=mycolor[2], xlab="H3K4me3 breadth (kb)",
        main="Distribution of H3K4me3 ChIP-seq peaks in NE2",border = NA,  font.lab = 2, cex.lab = 1.2)
dev.off()

#H3K4me1
NE2_K4me1 <- read.table("./data/ChIP/NE2H3K4me1.length.peak.bed",sep = "\t",header = F)
NE2_K4me1_2 <- NE2_K4me1[which(NE2_K4me1$V4 <= 25 ),]
x <- NE2_K4me1_2$V4
mycolor<-brewer.pal(9,"Set1")
pdf(file="Distribution of H3K4me1 ChIP-seq peaks in NE21.pdf")
h<-hist(x, breaks = 100, col=mycolor[2], xlab="H3K4me1 breadth (kb)",
        main="Distribution of H3K4me1 ChIP-seq peaks in NE2",border = NA,  font.lab = 2, cex.lab = 1.2)
dev.off()




