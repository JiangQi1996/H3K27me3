library(ggplot2)
setwd("D:/plot/figure3/TAD")
# matrix<- read.table("all.txt",sep = "\t",header = T)
# matrix$name <- factor(matrix$name,levels = c("broad","shuffle","narrow","control"))

broad_m <- read.table("broad_m.bed",sep = "\t",header = T)
# shuffle_m <- read.table("shuffle_m.txt",sep = "\t",header = T)
narrow_m <- read.table("narrow_m.bed",sep = "\t",header = T)
control_m <- read.table("control_m.bed",sep = "\t",header = T)

broad_m$number <- (broad_m$number)/1206  #落在TAD上的peak数目
narrow_m$number <- (narrow_m$number)/1182
control_m$number <- (control_m$number)/1291

library(RColorBrewer)
mycolor<- brewer.pal(9,"Set1")
colors <- c('darkred','firebrick','darkblue',"darkgrey")

library(ggplot2)
p1 <- ggplot(data = broad_m, mapping = aes(x = factor(group), y = number )) + geom_bar(stat = 'identity', position = 'dodge' , color = 'darkred' ,fill= "white",size = 1 )
# p2 <- p1 + geom_bar(data = shuffle_m, mapping = aes(x = factor(group), y = count,fill = 'firebrick'),stat = 'identity', position = 'dodge' , fill = 'darkblue')
p2 <- p1 + geom_bar(data = control_m, mapping = aes(x = factor(group), y = number,fill = 'darkblue'),stat = 'identity', position = 'dodge', color = 'darkgrey',fill="white", size = 1 )
p3 <- p2 + geom_bar(data = narrow_m, mapping = aes(x = factor(group), y = number,fill = 'darkgrey'),stat = 'identity', position = 'dodge', color = 'darkblue',fill= "white" ,size = 1 )
p4 = p3 + xlab('') + ylab('Number of H3K27me3 broad peaks') 
p4
#字体，坐标轴设置
p5 = p4 + theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), 
                panel.background = element_blank(),
                strip.background = element_blank(),
                axis.line = element_line(colour = "black",size = 1),
                title = element_text(face = "bold", colour = "black", size = 12),
                axis.title.x = element_text(size=12,face="bold",color="black",hjust=0.5),
                axis.title.y = element_text(size=12,face="bold",color="black",hjust=0.5),
                axis.text.x  = element_blank(),
                axis.text.y  = element_text(size = 9, color = "black", face = "bold"))
p5
ggsave(p5, file="TAD_hist.pdf", width=12, height=8)
