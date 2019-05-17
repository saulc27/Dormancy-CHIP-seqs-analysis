library(DESeq2)
library(ggplot2)
library(pheatmap)
library(pathfindR)
library(data.table)
library(dplyr)
library(fgsea)
library(tidyverse)
library(clusterProfiler)
library(plyr)
library(magrittr)
library(RColorBrewer)
library(DESeq)
library(gghighlight)
library(superheat)
library(ggpubr)
library(glimma)
library(ggiraph)
library(shiny)


setwd(dir = "../../Documents/School/BernsteinLab/Genomics/Chipseq/Dhep Thep/beds/Analysis")
read.promoter.counts <- read.table(file = "Promoters-with-DHEP3-THEP3-H3K27ac-and-BM-NR2F1-counts.csv", sep = ",", header = TRUE)
read.rnaseq.counts <- read.table(file = "DvT.RNAseq.results.csv", sep = ",", header = TRUE)
read.all.enhacers.table <- read.table(file = "../DHEP-THEP-H3K27ac-merged-q1e-3_filtered_AllEnhancers.table.txt", header = TRUE)
read.enhancer.promoter.table <- read.table(file ="../Dhep-Thep-enhancer-to-promoter-5-closest.bed", header = F)
all.enhancers.table.final <- read.all.enhacers.table[,c(1,2,3,4,10,12,13,14)]

CHIPseq.RNASEQ.counts.final <- merge(read.promoter.counts, read.rnaseq.counts, by.x="id", by.y="rn")
rownames(CHIPseq.RNASEQ.counts.final) = make.names(CHIPseq.RNASEQ.counts.final$id, unique=TRUE)

#Expression CHIP signal correlation
heatmap.plot.table <- CHIPseq.RNASEQ.counts.final %>% select(DHEP.K27ac, THEP.K27ac, NR2F1, baseMean)
t.heatmap.plot.table <- t(heatmap.plot.table)
superheat(t.heatmap.plot.table[c(1,2,3),],
          scale = F,
          yt = log2(t.heatmap.plot.table[4,]), yt.plot.type = "line",
          order.cols = order(t.heatmap.plot.table[4,], decreasing = T),
          #heat.col.scheme = "red",
          yt.axis.name = "baseMean counts")

#scatter plot log2fc vs chipsignal at promoters
scatter.plot.table <- CHIPseq.RNASEQ.counts.final %>% select(DHEP.K27ac, THEP.K27ac, NR2F1, log2FoldChange)
scatter.plot.table %>%
        mutate(highlight_flag = ifelse(abs(scatter.plot.table$log2FoldChange) > 2, T, F)) %>%
        ggplot(aes(x=log2(scatter.plot.table$NR2F1+2), y=scatter.plot.table$log2FoldChange)) +
        geom_point(aes(color = highlight_flag)) +
        scale_color_manual(values = c('#595959', 'red'))+
        #geom_point(size=1, color="#191970")+
        #geom_smooth(method=lm, color="#708090")+
        theme_minimal()+
        #scale_color_brewer(palette="Dark2")+
        xlab("Log2-Counts-NR2F1")+
        ylab("Log2FC")

#scatter plot H3K27ac DHEp3 v THEP3 at promoters
scatter.plot.table %>%
  mutate(highlight_flag = ifelse(abs(scatter.plot.table$log2FoldChange) > 2, T, F)) %>%
  ggplot(aes(x=log2(scatter.plot.table$THEP.K27ac+1), y=log2(scatter.plot.table$DHEP.K27ac+1))) +
  #geom_smooth(color="#8b5058")+ 
  #scale_color_manual(values = c('#595959', 'red'))+
  geom_point(aes(color = highlight_flag))+
  #geom_smooth(method=lm, color="#708090")+
  theme_minimal()+
  scale_color_manual(values = c('#595959', 'red'))+
  #scale_color_brewer(palette="Blues")+
  xlab("Log2-Counts-H3K27ac-THEP3")+
  ylab("Log2-Counts-H3K27ac-DHEP3")

#interactive
setDT(scatter.plot.table, keep.rownames = TRUE)[]
my_gg <- g + geom_point_interactive(aes(tooltip = rn), size = 2) 
#girafe(code = print(my_gg) )







#boxplot
log2FC <- scatter.plot.table %>% filter(log2FoldChange>2)
log2FC$category<-"UP"
log2FC <- log2FC %>%  select(DHEP.K27ac, THEP.K27ac, NR2F1, category)
up.table <- melt(log2FC, id.var = "category")

logneg2FC <- scatter.plot.table %>% filter(log2FoldChange< -2) 
logneg2FC$category<-"Down"
logneg2FC <- logneg2FC %>%  select(DHEP.K27ac, THEP.K27ac, NR2F1, category)
down.table <- melt(logneg2FC, id.var = "category")

up.down.table <- rbind(up.table, down.table) 

ggplot(up.down.table,aes(x=variable,y=log2(value), fill=category))+ #facet_wrap(~variable)+
  theme_bw() +  scale_fill_brewer(palette="Blues")+
  geom_boxplot(notch=TRUE, alpha=1, width=0.5, size=0, position=position_dodge(0.65))+ 
  stat_compare_means(label =  "p.signif", label.x = 1.5, method = "t.test")
  #geom_violin(trim = FALSE, add = "boxplot")
    
ggdensity(up.down.table, x = "value",
          add = "mean", rug = TRUE,
          color = "variable", fill = "variable",
          palette = c("#00AFBB", "#E7B800", "#708090"))+ xscale("log2") 


