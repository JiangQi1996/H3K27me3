#计算multiple-interactioon 的peak number
setwd("D:/plot/figure3/mutiplex-chromatin-interaction")
broad_m <- read.table("broad_anchor_num.txt",sep="\t",header = T)
broad_unique <- unique(broad_m[,c(4,9)])
broad_fre <- table(broad_unique$peak_id)
broad_fre <- as.data.frame(broad_fre)
nrow(broad_fre[which(broad_fre$Freq >1),])
levels(factor(broad_fre[,2]))
write.table(broad_fre,file ="broad_fre.txt", sep = "\t",row.names = FALSE, col.names = TRUE, quote =FALSE)

narrow_m <- read.table("narrow_anchor_num.txt",sep="\t",header = T)
narrow_unique <- unique(narrow_m[,c(4,9)])
narrow_fre <- table(narrow_unique$peak_id)
narrow_fre <- as.data.frame(narrow_fre)
nrow(narrow_fre[which(narrow_fre$Freq >1),])
levels(factor(narrow_fre[,2]))
write.table(narrow_fre,file ="narrow_fre.txt", sep = "\t",row.names = FALSE, col.names = TRUE, quote =FALSE)

control_m <- read.table("control_anchor_num.txt",sep="\t",header = T)
control_unique <- unique(control_m[,c(4,9)])
control_fre <- table(control_unique$peak_id)
control_fre <- as.data.frame(control_fre)
nrow(control_fre[which(control_fre$Freq >1),])
levels(factor(control_fre[,2]))
write.table(control_fre,file ="control_fre.txt", sep = "\t",row.names = FALSE, col.names = TRUE, quote =FALSE)

#统计broad peak落在不同数目的anchor的比例得到interaction_num.txt文件
nrow(broad_fre[which(broad_fre$Freq == i),]) / length(broad_fre$Var1)

#R plot 
library(RColorBrewer)
library(ggplot2)
setwd("D:/plot/figure3/mutiplex-chromatin-interaction")
interaction_num <- read.table("interaction_num.txt",sep="\t",header = T,row.names= 1)
colors <- c('darkred','navyblue',"darkgrey")
interaction_num <- as.matrix(interaction_num)
interaction_num <- t(interaction_num)
pdf("interaction_frequency_barplot.pdf", width = 8, height = 8)
a <-barplot(
  interaction_num,
  col = colors,
  xlab = expression(bold('Interaction number')),
  ylab = expression(bold('H3K27me3 peaks')),
  border = 'black',
  cex.lab = 1.5,
  font.lab = 2,
  axes = F,
  width = c(0.5,0.5,0.5,0.5),
  names.arg = c("1","2","3",">=4"),
  beside=TRUE,
  ylim = c(0,0.7)
)
axis(2,lwd=3,lwd.ticks = 3,cex.axis=1.5,font.axis=2,col="black")
# legend(3.2,280, rownames(Values), cex = 1.3, fill = colors)
# text(aes(label=paste(Values)),angle=15,vjust=-0.2,position=position_dodge(0.9),size=5)
# text(x = a, y = Values, label=Values, cex=1, pos=3, col="black", xpd=TRUE)
dev.off()