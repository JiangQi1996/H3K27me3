#计算peak之间interaction number
setwd("D:/plot/figure3/peak-peak-interaction")
length <- read.table("loop_length.bed",sep="\t",header=F)
colnames(length) <- c("loop_id","length")
broad_m <- read.table("broad_anchor.bed",sep="\t",header=F)
broad_m[,10] <- paste(paste(paste(broad_m[,1],broad_m[,2],sep="-"),broad_m[,3],sep="-"),broad_m[,4],sep="-")
B <- unique(broad_m[,1:4])
B[,5] = rep("1",nrow(B))
B[,6] <- paste(paste(B[,1],B[,2],sep="-"),B[,3],sep="-")
B[,5] <- as.numeric(B[,5])
C <- aggregate(B[,5], by = list(B[,6]), FUN=sum)
library(plyr)
colnames(C) = c("Class", "Num")
colnames(B) = c("S1","S2","S3","S4","num","Class")
D <- join(B, C,by = "Class", type = "left")
D[,8] <- paste(D[,6],D[,4],sep="-")
colnames(D) = c("S1","S2","S3","S4","num","Class","Num","Class2")
colnames(broad_m) = c("S1","S2","S3","loop_id","S5","S6","S7","S8","S9","Class2")
E <- join(broad_m, D,by = "Class2", type = "left")
E <- E[,c(1:9,17)]
E<- merge(E,length,by="loop_id")
E <- E[,c(2,3,4,1,11,5,6,7,8,9,10)]
colnames(E) <- c("chr1","start1","end1","loop_id","loop_length(Mb)","chr2","start2","end2","peak_id","peak_length(kb)","interaction_number")
write.table(E,file ="broad_anchor_num.txt", sep = "\t",row.names = FALSE, col.names = TRUE, quote =FALSE)

narrow_m <- read.table("narrow_anchor.bed",sep="\t",header=F)
narrow_m[,10] <- paste(paste(paste(narrow_m[,1],narrow_m[,2],sep="-"),narrow_m[,3],sep="-"),narrow_m[,4],sep="-")
B <- unique(narrow_m[,1:4])
B[,5]=rep("1",nrow(B))
B[,6] <- paste(paste(B[,1],B[,2],sep="-"),B[,3],sep="-")
B[,5] <- as.numeric(B[,5])
C <- aggregate(B[,5], by = list(B[,6]), FUN=sum)
library(plyr)
colnames(C) = c("Class", "Num")
colnames(B) = c("S1","S2","S3","S4","num","Class")
D <- join(B, C,by = "Class", type = "left")
D[,8] <- paste(D[,6],D[,4],sep="-")
colnames(D) = c("S1","S2","S3","S4","num","Class","Num","Class2")
colnames(narrow_m) = c("S1","S2","S3","loop_id","S5","S6","S7","S8","S9","Class2")
E <- join(narrow_m, D,by = "Class2", type = "left")
E <- E[,c(1:9,17)]
E<- merge(E,length,by="loop_id")
E <- E[,c(2,3,4,1,11,5,6,7,8,9,10)]
colnames(E) <- c("chr1","start1","end1","loop_id","loop_length(Mb)","chr2","start2","end2","peak_id","peak_length(kb)","interaction_number")
write.table(E,file ="narrow_anchor_num.txt", sep = "\t",row.names = FALSE, col.names = TRUE, quote =FALSE)

control_m <- read.table("control_anchor.bed",sep="\t",header=F)
control_m[,10] <- paste(paste(paste(control_m[,1],control_m[,2],sep="-"),control_m[,3],sep="-"),control_m[,4],sep="-")
B <- unique(control_m[,1:4])
B[,5]=rep("1",nrow(B))
B[,6] <- paste(paste(B[,1],B[,2],sep="-"),B[,3],sep="-")
B[,5] <- as.numeric(B[,5])
C <- aggregate(B[,5], by = list(B[,6]), FUN=sum)
library(plyr)
colnames(C) = c("Class", "Num")
colnames(B) = c("S1","S2","S3","S4","num","Class")
D <- join(B, C,by = "Class", type = "left")
D[,8] <- paste(D[,6],D[,4],sep="-")
colnames(D) = c("S1","S2","S3","S4","num","Class","Num","Class2")
colnames(control_m) = c("S1","S2","S3","loop_id","S5","S6","S7","S8","S9","Class2")
E <- join(control_m, D,by = "Class2", type = "left")
E <- E[,c(1:9,17)]
E<- merge(E,length,by="loop_id")
E <- E[,c(2,3,4,1,11,5,6,7,8,9,10)]
colnames(E) <- c("chr1","start1","end1","loop_id","loop_length(Mb)","chr2","start2","end2","peak_id","peak_length(kb)","interaction_number")
write.table(E,file ="control_anchor_num.txt", sep = "\t",row.names = FALSE, col.names = TRUE, quote =FALSE)

#R plot for peak-peak-interaction
setwd("D:/plot/figure3/peak-peak-interaction")
broad_m <- read.table("broad_anchor_num.txt",sep="\t",header = T)
broad_unique <- unique(broad_m[,c(4,9)])
broadloop_fre <- table(broad_unique$loop_id)
broadloop_fre <- as.data.frame(broadloop_fre)
a <- broadloop_fre[which(broadloop_fre$Freq >1),1]
a<-as.data.frame(a)
colnames(a) <- c("loop_id")
broad_loop_peak_num <- nrow(merge(a,broad_unique))
broad_peak_num <- unique(broad_unique$peak_id)

narrow_m <- read.table("narrow_anchor_num.txt",sep="\t",header = T)
narrow_unique <- unique(narrow_m[,c(4,9)])
narrowloop_fre <- table(narrow_unique$loop_id)
narrowloop_fre <- as.data.frame(narrowloop_fre)
a <- narrowloop_fre[which(narrowloop_fre$Freq >1),1]
a<-as.data.frame(a)
colnames(a) <- c("loop_id")
narrow_loop_peak_num <- nrow(merge(a,narrow_unique))
narrow_peak_num <- unique(narrow_unique$peak_id)

control_m <- read.table("control_anchor_num.txt",sep="\t",header = T)
control_unique <- unique(control_m[,c(4,9)])
controlloop_fre <- table(control_unique$loop_id)
controlloop_fre <- as.data.frame(controlloop_fre)
a <- controlloop_fre[which(controlloop_fre$Freq >1),1]
a<-as.data.frame(a)
colnames(a) <- c("loop_id")
control_loop_peak_num <- nrow(merge(a,control_unique))
control_peak_num <- unique(control_unique$peak_id)

value <- c(broad_loop_peak_num/431,2/47,4/84)
pdf("peak proportion between loops.pdf", width = 4, height = 4)
mycolor <- c('darkred','darkblue',"darkgrey")
barplot(
  value,
  col = mycolor,
  ylab = expression(bold('Peak Proportion of interaction')),
  border = mycolor,
  cex.lab = 1.5,
  font.lab = 2,
  axes = F,
  width = c(0.2,0.2,0.2),
  space = 0.5,
  names.arg = c("broad","narrow","control")
)
axis(2,lwd=3,lwd.ticks = 3,cex.axis=1.5,font.axis=2,col="black")
# legend(3.2,800, c("broad","narrow","random"),fill = c(mycolor[c(1,2,3)]))
dev.off()

