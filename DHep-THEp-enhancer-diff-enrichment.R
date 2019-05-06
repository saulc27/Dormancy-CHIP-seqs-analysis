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
library(ggstatsplot)
library(superheat)

setwd(dir = "../../Documents/School/BernsteinLab/Genomics/Chipseq/Dhep Thep/beds/")
read.k27ac.rose.peaks <- read.table(file = "DHEP-THEP-H3K27ac-rose-peaks-table.csv", sep = ",", header = TRUE)
colnames(read.k27ac.rose.peaks) <- c("Peak-ID","THEP3-H3K27ac","DHEP3-H3K27ac")
rose.peaks.table <- read.k27ac.rose.peaks %>%
                    `row.names<-`(., NULL) %>% 
                     column_to_rownames(var = "Peak-ID")
                

#rose.peaks.table.ratio <- k27ac.rose.peaks %>%
 # `row.names<-`(., NULL) %>% 
  #column_to_rownames(var = "Peak-ID")%>%
  #mutate(THEP3-H3K27ac = THEP3-H3K27ac+1)%>%
  #mutate(DHEP3-H3K27ac = DHEP3-H3K27ac+1)
  #mutate(ratio = `DHEP3-H3K27ac`/`THEP3-H3K27ac`)


rose.peaks.table$`THEP3-H3K27ac` <- as.integer(rose.peaks.table$`THEP3-H3K27ac`)
rose.peaks.table$`DHEP3-H3K27ac` <- as.integer(rose.peaks.table$`DHEP3-H3K27ac`)
            
condition <- factor(c("THEP3","DHEP3"))
coldata <- data.frame(row.names=colnames(rose.peaks.table), condition)
#dds <- DESeqDataSetFromMatrix(countData=rose.peaks.table, colData=coldata, design=1)
dds <- DESeqDataSetFromMatrix(countData=rose.peaks.table, colData=coldata, design=~condition)

dds$condition <- relevel(dds$condition, "DHEP3")
dds <- estimateSizeFactors(dds)
vst <- vst(dds, blind = TRUE)
ggplot(rose.peaks.table, aes(x=log2((rose.peaks.table$`THEP3-H3K27ac`)+1), y=log2((rose.peaks.table$`DHEP3-H3K27ac`)+1))) +
                            geom_point(size=0.7, color="#191970")+
                            geom_smooth(method=lm, color="#708090")+
                            theme_minimal()+
                            scale_color_brewer(palette="Dark2")+
                            xlab("THEP3-H3K27ac")+
                            ylab("DHEP3-H3K27ac")

dds <- DESeq(dds)
res<- results(dds)
table(res$padj<0.05)
res.table <- data.frame(res)

promoter.merge.table <- read.csv(file = "../../../RNAseq R/Dormancy/RNAseq-promoter-merge-with-counts.csv")

BM.table <- read.csv(file = "../../../RNAseq R/Dormancy/Norm.counts.BM.DHEP.THEP.csv", row.names = "X")
keep <- rowSums(BM.table) > 400
BM.table.keep <- BM.table[keep,]

promoter.signal.NR2F1. <- read.csv(file = "../../../RNAseq R/Dormancy/Promoters-with-DHEP3-THEP3-H3K27ac-and-BM-NR2F1-counts.csv")

 


merge.symbol.signal.fc <- data.frame(promoter.merge.table$Symbol, data.frame(promoter.merge.table$X...LOG10qvalue.NR2F1, data.frame(promoter.merge.table$log2FoldChange)))
colnames(merge.symbol.signal.fc) <- c("Symbol", "qVal", "Log2FC")

#add color to selected FC points                                     
merge.symbol.signal.fc %>%
  mutate(highlight_flag = ifelse(log2(merge.symbol.signal.fc$qVal) > 2  & abs(merge.symbol.signal.fc$Log2FC) > 2, T, F)) %>%
      ggplot(aes(x=log2(merge.symbol.signal.fc$qVal), y=merge.symbol.signal.fc$Log2FC)) +
        geom_point(aes(color = highlight_flag)) +
        scale_color_manual(values = c('#595959', 'red'))+
      #geom_point(size=1, color="#191970")+
      #geom_smooth(method=lm, color="#708090")+
        theme_minimal()+
      #scale_color_brewer(palette="Dark2")+
        xlab("Log2-qVal-NR2F1")+
        ylab("Log2FC")
     
#ggstatsplot DHEP3 replicates
#tiff("DHEP-replicates-scatter-plot-ggstatsplot.tiff", units="in", width=5, height=5, res=300)
ggscatterstats(
  data = log2(rose.peaks.table+1),
  x = 1,
  y = 2,
  xlab = "Log2 Counts THEP3 H3K27ac",
  ylab = "Log2 Counts DHEP3 H3K27ac",
  title = "H3K27ac CHIP-seq counts between DHEP3 and THEP3",
  point.size = 1, point.alpha = 1, line.size =0.5,
  bf.message = FALSE,
  messages = FALSE,
  marginal = TRUE,
  method = "lm",
  results.subtitle = TRUE
)

BM.exp.NR2F1.merge <- merge(promoter.signal.NR2F1., BM.table.keep, by.x="id", by.y="row.names")

mat.TNF.LE <- merge(assay(rld), Hallmark.TNF.LE, by.x="row.names", by.y="genes")

promoter.merge.table.sorted <- promoter.merge.table %>%
  arrange(desc(baseMean))

#promoter.merge.table.sorted$baseMean <- as.character(promoter.merge.table.sorted$baseMean)
#promoter.merge.table.sorted$baseMean  <- factor(promoter.merge.table.sorted$baseMean , levels=unique(promoter.merge.table.sorted$baseMean ))


#mean vs signal plot
heatmap.plot.table <- promoter.merge.table %>% select(Symbol, X...LOG10qvalue.NR2F1, baseMean, padj)
#make rownames unique
rownames(heatmap.plot.table) = make.names(heatmap.plot.table$Symbol, unique=TRUE)

heatmap.plot.table <-heatmap.plot.table %>% select(-Symbol)
colnames(heatmap.plot.table) <- c("qValue","baseMean","padj")
t.heatmap.plot.table <- t(heatmap.plot.table)

heatmap.plot.BM <- BM.NR2F1.merge %>% select(id, DHEP.K27ac, THEP.K27ac, NR2F1, BM3, BM4) 
rownames(heatmap.plot.BM) = make.names(heatmap.plot.BM$id, unique=TRUE)
heatmap.plot.BM <-heatmap.plot.BM %>% select(-id)
heatmap.plot.BM.mat <- as.matrix(heatmap.plot.BM)
t.heatmap.plot.table.BM <- t(heatmap.plot.BM.mat)


heatmap.plot.BM.mat.log <-  heatmap.plot.BM.mat+1
t.heatmap.plot.table.BM <- t(heatmap.plot.BM.mat.log)


superheat(t.heatmap.plot.table.BM[c(3,4),],
          scale = F,
          yt = log2(t.heatmap.plot.table.BM[4,]), yt.plot.type = "line",
          order.cols = order(t.heatmap.plot.table.BM[4,], decreasing = T),
          yt.axis.name = "baseMean counts")


         
          