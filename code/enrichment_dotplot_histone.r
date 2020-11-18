##enrichment dotplot of NE2 histone

library(ggplot2)
library(RColorBrewer)
mycolor<-brewer.pal(9,"Blues")

#H3K27me3
pathway <- read.table("./data/ChIP/H3K27me3.txt",sep="\t",header = T)
H3K27me3_raw <- read.csv("./data/ChIP/GO_NE2_H3K27me3.csv",header = T)
H3K27me3_pathway <- merge(pathway,H3K27me3_raw,by="Description",sort = F)
write.table(H3K27me3_pathway,file ="./result/H3K27me3_pathway.txt", sep = "\t",row.names = FALSE, col.names = TRUE, quote =FALSE)
H3K27me3_pathway$LogP <- abs(H3K27me3_pathway$LogP)
a<-unique(H3K27me3_pathway$Description)
a <- rev(a)
H3K27me3_pathway$Description <- factor(H3K27me3_pathway$Description,levels = a)
H3K27me3_pathway$LogP[H3K27me3_pathway$LogP > 10]<-10
p1 = ggplot(H3K27me3_pathway,aes(GeneList,Description,color = LogP,size = Enrichment)) + geom_point(shape = 20)
p1 = p1 + scale_colour_gradient(high =  mycolor[8],low = "white") + xlab("H3K27me3 breadth") + ylab("")
p1 = p1 + labs(color=expression(-log[10](Pvalue)),size="Fold enrichment")
p1 = p1 + theme(panel.background = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line = element_line(colour = "black",size = 0.8),
                axis.text.x = element_blank(),
                axis.text.y = element_text(size = 12, color = "black", face = "plain"),
                axis.title.x = element_text(size = 12, color = "black", face = "bold"),
                axis.ticks = element_blank(),
                panel.border = element_rect(colour="black",fill=NA,size=1))
p1 = p1 + scale_size(limits=c(1,12),range = c(-1,10))
ggsave(plot = p4, height=8, width=10, dpi=600, filename="./result/H3K27me3.dotplot.pdf", useDingbats=FALSE)

#H3K27ac
pathway <- read.table("./data/ChIP/H3K27ac.txt",sep="\t",header = T)
H3K27ac_raw <- read.csv("./data/ChIP/GO_NE2_H3K27ac.csv",header = T)
H3K27ac_pathway <- merge(pathway,H3K27ac_raw,by="Description",sort = F)
write.table(H3K27ac_pathway,file ="./result/H3K27ac_pathway.txt", sep = "\t",row.names = FALSE, col.names = TRUE, quote =FALSE)
H3K27ac_pathway$LogP <- abs(H3K27ac_pathway$LogP)
a<-unique(H3K27ac_pathway$Description)
a <- rev(a)
H3K27ac_pathway$Description <- factor(H3K27ac_pathway$Description,levels = a)
H3K27ac_pathway$LogP[H3K27ac_pathway$LogP > 10]<-10
p2 = ggplot(H3K27ac_pathway,aes(GeneList,Description,color = LogP,size = Enrichment)) + geom_point(shape = 20)
p2 = p2 + scale_colour_gradient(high =  mycolor[8],low = "white") + xlab("H3K27ac breadth") + ylab("")
p2 = p2 + labs(color=expression(-log[10](Pvalue)),size="Fold enrichment")
p2 = p2 + theme(panel.background = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line = element_line(colour = "black",size = 0.8),
                axis.text.x = element_blank(),
                axis.text.y = element_text(size = 12, color = "black", face = "plain"),
                axis.title.x = element_text(size = 12, color = "black", face = "bold"),
                axis.ticks = element_blank(),
                panel.border = element_rect(colour="black",fill=NA,size=1))
p2 = p2 + scale_size(limits=c(1,4),range = c(-1,8))
ggsave(plot = p2, height=8, width=10, dpi=600, filename="./result/H3K27ac.dotplot.pdf", useDingbats=FALSE)

#H3K4me3
pathway <- read.table("./data/ChIP/H3K4me3.txt",sep="\t",header = T)
H3K4me3_raw <- read.csv("./data/ChIP/GO_NE2_H3K4me3.csv",header = T)
H3K4me3_pathway <- merge(pathway,H3K4me3_raw,by="Description",sort = F)
write.table(H3K4me3_pathway,file ="./result/H3K4me3_pathway.txt", sep = "\t",row.names = FALSE, col.names = TRUE, quote =FALSE)
H3K4me3_pathway$LogP <- abs(H3K4me3_pathway$LogP)
a<-unique(H3K4me3_pathway$Description)
a <- rev(a)
H3K4me3_pathway$Description <- factor(H3K4me3_pathway$Description,levels = a)
H3K4me3_pathway$LogP[H3K4me3_pathway$LogP > 10]<-10
p3 = ggplot(H3K4me3_pathway,aes(GeneList,Description,color = LogP,size = Enrichment)) + geom_point(shape = 20)
p3 = p3 + scale_colour_gradient(high = mycolor[8],low = "white") + xlab("H3K4me3 breadth") + ylab("")
p3 = p3 + labs(color=expression(-log[10](Pvalue)),size="Fold enrichment")
p3 = p3 + theme(panel.background = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line = element_line(colour = "black",size = 0.8),
                axis.text.x = element_blank(),
                axis.text.y = element_text(size = 12, color = "black", face = "plain"),
                axis.title.x = element_text(size = 12, color = "black", face = "bold"),
                axis.ticks = element_blank(),
                panel.border = element_rect(colour="black",fill=NA,size=1))
p3 = p3 + scale_size(limits=c(1,6),range = c(-1,8))
ggsave(plot = p3, height=8, width=10, dpi=600, filename="./result/H3K4me3.dotplot.pdf", useDingbats=FALSE)

#H3K4me1
pathway <- read.table("./data/ChIP/H3K4me1.txt",sep="\t",header = T)
H3K4me1_raw <- read.csv("./data/ChIP/GO_NE2_H3K4me1.csv",header = T)
H3K4me1_pathway <- merge(pathway,H3K4me1_raw,by="Description",sort = F)
write.table(H3K4me1_pathway,file ="./result/H3K4me1_pathway.txt", sep = "\t",row.names = FALSE, col.names = TRUE, quote =FALSE)
H3K4me1_pathway$LogP <- abs(H3K4me1_pathway$LogP)
a<-unique(H3K4me1_pathway$Description)
a <- rev(a)
H3K4me1_pathway$Description <- factor(H3K4me1_pathway$Description,levels = a)
H3K4me1_pathway$LogP[H3K4me1_pathway$LogP > 10]<-10
p4 = ggplot(H3K4me1_pathway,aes(GeneList,Description,color = LogP,size = Enrichment)) + geom_point(shape = 20)
p4 = p4 + scale_colour_gradient(high =  mycolor[8],low = "white") + xlab("H3K4me1 breadth") + ylab("")
p4 = p4 + labs(color=expression(-log[10](Pvalue)),size="Fold enrichment")
p4 = p4 + theme(panel.background = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line = element_line(colour = "black",size = 0.8),
                axis.text.x = element_blank(),
                axis.text.y = element_text(size = 12, color = "black", face = "plain"),
                axis.title.x = element_text(size = 12, color = "black", face = "bold"),
                axis.ticks = element_blank(),
                panel.border = element_rect(colour="black",fill=NA,size=1))
p4 = p4 + scale_size(limits=c(1,5),range = c(-1,8))
p4
ggsave(plot = p4, height=8, width=10, dpi=600, filename="./result/H3K4me1.dotplot.pdf", useDingbats=FALSE)


#IMR90
pathway <- read.table("./data/ChIP/IMR90.txt",sep="\t",header = T)
H3K27me3_raw <- read.csv("./data/ChIP/GO_IMR90_H3K27me3.csv",header = T)
H3K27me3_pathway <- merge(pathway,H3K27me3_raw,by="Description",sort = F)
write.table(H3K27me3_pathway,file ="./result/IMR90_H3K27me3_pathway.txt", sep = "\t",row.names = FALSE, col.names = TRUE, quote =FALSE)
H3K27me3_pathway$LogP <- abs(H3K27me3_pathway$LogP)
a<-unique(H3K27me3_pathway$Description)
a <- rev(a)
H3K27me3_pathway$Description <- factor(H3K27me3_pathway$Description,levels = a)
H3K27me3_pathway$LogP[H3K27me3_pathway$LogP > 10]<-10
p5 = ggplot(H3K27me3_pathway,aes(GeneList,Description,color = LogP,size = Enrichment)) + geom_point(shape = 20)
p5 = p5 + scale_colour_gradient(high =  mycolor[8],low = "white") + xlab("H3K27me3 breadth") + ylab("")
p5 = p5 + labs(color=expression(-log[10](Pvalue)),size="Fold enrichment")
p5 = p5 + theme(panel.background = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line = element_line(colour = "black",size = 0.8),
                axis.text.x = element_blank(),
                axis.text.y = element_text(size = 12, color = "black", face = "plain"),
                axis.title.x = element_text(size = 12, color = "black", face = "bold"),
                axis.ticks = element_blank(),
                panel.border = element_rect(colour="black",fill=NA,size=1))
p5 = p5 + scale_size(limits=c(1,11),range = c(1,10))
ggsave(plot = p5, height=8, width=10, dpi=600, filename="./result/IMR90.H3K27me3.dotplot.pdf", useDingbats=FALSE)

