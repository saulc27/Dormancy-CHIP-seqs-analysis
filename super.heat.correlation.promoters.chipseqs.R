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
#active promoters
read.promoter.counts <- read.table(file = "../Promoters-with-DHEP3-THEP3-H3K27ac-and-BM-NR2F1-counts.csv", sep = ",", header = TRUE)

#RNAseq table
read.rnaseq.results <- read.table(file = "../DvT.RNAseq.results.csv", sep = ",", header = TRUE)

#ROSe enhancer table
all.enhancers.table <- read.table(file = "../DHEP-THEP-H3K27ac-merged-q1e-3_filtered_AllEnhancers.table.txt", header = TRUE)[,c(1,2,3,4,10,12,13,14)]

#Enhancer promoter tables (bedtools)
read.1mb.enhancer.promoter.table <- read.table(file ="DHEP-THEP-H3K27ac-merged-q1e-3_filtered_Enhancers_withSuper-1mb-extended-with-promoters.bed", header = F)
read.2mb.enhancer.promoter.table <- read.table(file ="DHEP-THEP-H3K27ac-merged-q1e-3_filtered_Enhancers_withSuper-2mb-extended-with-promoters.bed", header = F)

#SE only with Log2Fc cutoff
SE.only.table <- all.enhancers.table[all.enhancers.table$isSuper==1,]
SE.only.table$log2FC  <- log2(SE.only.table$THEP_H3K27ac.XL.final.bam/SE.only.table$DHEP_H3K27ac.XL.final.bam)
SE.only.table.DF <- SE.only.table[which(SE.only.table$log2FC>1.5 | SE.only.table$log2FC< -1.5),]

#Merging DF SE enhancer table with promoters removing dups
SE.only.with.promoters <- merge(read.1mb.enhancer.promoter.table, SE.only.table.DF, by.x = "V4", by.y="REGION_ID")
SE.only.with.promoters <- SE.only.with.promoters[!duplicated(SE.only.with.promoters),]
Trimmed.SE.only.with.promoters <- SE.only.with.promoters[,c(-2,-3,-4,-5,-6,-7,-8)]

#merging with RNAseq table and filtering lower expressed genes
RNA.seq.SE.table.promoters<- merge(Trimmed.SE.only.with.promoters , read.rnaseq.results, by.x="V9", by.y="rn")
RNA.seq.SE.table.promoters.filtered<- RNA.seq.SE.table.promoters[RNA.seq.SE.table.promoters$baseMean>500,]

SE.table.promoters.RNA.seq.THEP <- RNA.seq.SE.table.promoters.filtered[RNA.seq.SE.table.promoters.filtered$log2FC>1,] 
SE.table.promoters.RNA.seq.THEP <- SE.table.promoters.RNA.seq.THEP[!duplicated(SE.table.promoters.RNA.seq.THEP),]

SE.table.promoters.RNA.seq.DHEP <-RNA.seq.SE.table.promoters.filtered[RNA.seq.SE.table.promoters.filtered$log2FC< -1,] 
SE.table.promoters.RNA.seq.DHEP <- SE.table.promoters.RNA.seq.DHEP[!duplicated(SE.table.promoters.RNA.seq.DHEP),]

Final.SE.table.promoters.RNA.seq.THEP <- SE.table.promoters.RNA.seq.THEP[SE.table.promoters.RNA.seq.THEP$log2FoldChange>0.8,]
write.csv(Final.SE.table.promoters.RNA.seq.THEP, file = "SE.genes.table.THEP3.csv")


final.table.promoters.RNA.seq.THEP.top <- SE.table.promoters.RNA.seq.THEP  %>%
  group_by(V4) %>%
  arrange(desc(log2FoldChange)) %>% 
  slice(1:20)


Final.SE.table.promoters.RNA.seq.DHEP <- SE.table.promoters.RNA.seq.DHEP[SE.table.promoters.RNA.seq.DHEP$log2FoldChange< -0.8,]
write.csv(Final.SE.table.promoters.RNA.seq.DHEP, file = "SE.genes.table.DHEP3.csv")


final.table.promoters.RNA.seq.DHEP.top <- SE.table.promoters.RNA.seq.DHEP  %>%
  group_by(V4) %>%
  arrange(log2FoldChange) %>% 
  slice(1:5)



#TE only with Log2FC
TE.only.table <- all.enhancers.table[all.enhancers.table$isSuper==0,]
TE.only.table$log2FC  <- log2((TE.only.table$THEP_H3K27ac.XL.final.bam+1)/(TE.only.table$DHEP_H3K27ac.XL.final.bam+1))
TE.only.table$CENTER <- as.integer((TE.only.table$START+TE.only.table$STOP)/2) 
TE.only.table.DF <- TE.only.table[which(TE.only.table$log2FC>1.5 | TE.only.table$log2FC< -1.5),]


#Merging DF TE enhancer table with promoters removing dups
TE.only.with.promoters <- merge(read.1mb.enhancer.promoter.table, TE.only.table.DF, by.x = "V4", by.y="REGION_ID")
TE.only.with.promoters <- TE.only.with.promoters[!duplicated(TE.only.with.promoters),]
Trimmed.TE.only.with.promoters <- TE.only.with.promoters[,c(-2,-3,-4,-5)]
Trimmed.TE.only.with.promoters <- Trimmed.TE.only.with.promoters[Trimmed.TE.only.with.promoters$V7>0,]
Trimmed.TE.only.with.promoters$TSS <- as.integer((Trimmed.TE.only.with.promoters$V7+Trimmed.TE.only.with.promoters$V8)/2)
Trimmed.TE.only.with.promoters$distance <- abs(Trimmed.TE.only.with.promoters$CENTER - Trimmed.TE.only.with.promoters$TSS)

#adding ranks by distance
Trimmed.TE.only.with.promoters <-   Trimmed.TE.only.with.promoters %>%
                                    group_by(V4) %>%
                                    mutate(my_ranks = order(order(distance, decreasing=F)))


#merging with RNAseq table and filtering lower expressed genes
RNA.seq.TE.table.promoters<- merge(Trimmed.TE.only.with.promoters , read.rnaseq.results, by.x="V9", by.y="rn")
RNA.seq.TE.table.promoters.filtered<- RNA.seq.TE.table.promoters[RNA.seq.TE.table.promoters$baseMean>250,]
TE.table.promoters.RNA.seq.THEP <- RNA.seq.TE.table.promoters.filtered[RNA.seq.TE.table.promoters.filtered$log2FC>1,] 
TE.table.promoters.RNA.seq.THEP <- TE.table.promoters.RNA.seq.THEP[!duplicated(TE.table.promoters.RNA.seq.THEP),]
TE.table.promoters.RNA.seq.DHEP <-RNA.seq.TE.table.promoters.filtered[RNA.seq.TE.table.promoters.filtered$log2FC< -1,] 
TE.table.promoters.RNA.seq.DHEP <- TE.table.promoters.RNA.seq.DHEP[!duplicated(TE.table.promoters.RNA.seq.DHEP),]

Final.TE.table.promoters.RNA.seq.DHEP <- TE.table.promoters.RNA.seq.THEP[TE.table.promoters.RNA.seq.THEP$log2FoldChange> 0.8,]
write.csv(Final.TE.table.promoters.RNA.seq.DHEP, file = "TE.genes.table.THEP3.csv")



TE.final.table.promoters.RNA.seq.THEP.top <- TE.table.promoters.RNA.seq.THEP  %>%
  group_by(V4) %>%
    arrange(log2FoldChange) %>% 
  slice(1)

TE.final.table.promoters.RNA.seq.DHEP.top <- TE.table.promoters.RNA.seq.DHEP  %>%
  group_by(V4) %>%
  arrange(log2FoldChange) %>% 
  slice(1)





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


