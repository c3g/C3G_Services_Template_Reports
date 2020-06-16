closeAllConnections()
rm(list = ls())

setwd(
  "/Users/skailasam/Desktop/mnt/projects/Petunia_Tobacco_RNAseq_C3G-IM-427/01_Petunia/Analysis_report/"
)

a <- c(1:100)
figcount <- rep(sprintf("Figure_%02d_", a))
cnt <- 0


# ##############################################################################################################################
figdir <-
  c(
    "/Users/skailasam/Desktop/mnt/projects/Petunia_Tobacco_RNAseq_C3G-IM-427/01_Petunia/Analysis_report/output_FC_2.0_FDR_0.05/"
  )
tabdir <-
  c(
    "/Users/skailasam/Desktop/mnt/projects/Petunia_Tobacco_RNAseq_C3G-IM-427/01_Petunia/Analysis_report/output_FC_2.0_FDR_0.05/Tables/"
  )
FCcutoff <- c(2)  # FDR < 5 %
log2FCcutval <- c(log2(FCcutoff))
aPVALcutoff <- c(0.05)  # adj.P.val < 5 %
avglogCPMcutoff <- 4
# ##############################################################################################################################


# myDir<- setwd("output/")
# file.remove(file.path(myDir, list.files()))



#library(limma)
library(edgeR)
#library(DESeq2)
# library(statmod)
# library(sva)
# library(qvalue)
library(RColorBrewer)
# library(plyr)
library(ggplot2)
library(reshape2)
library(gplots)
# library(Vennerable)
library(ggpubr)
library(gridExtra)
library(grid)
library(cowplot)
library(magick)
library(ComplexHeatmap)
library(data.table)

## ----readcounts----------------------------------------------------------
GenewiseCounts <- read.delim("/Users/skailasam/Desktop/mnt/projects/Petunia_Tobacco_RNAseq_C3G-IM-427/01_Petunia/DGE/rawCountMatrix.csv",
                             header = T,
                             sep = "\t",
                             row.names = 1
)
genes <- GenewiseCounts[c(1)]
colnames(GenewiseCounts) <-
  substring(colnames(GenewiseCounts), 1, 7)

#reorder
#GenewiseCounts<- GenewiseCounts[c("Symbol","DPP_CN1","DPP_CN2","DPP_CN3","DPP_CN4","DPP_CN5","DPP_EX1","DPP_EX2","DPP_EX3","DPP_EX4","DPP_EX5","PPP_CN1","PPP_CN2","PPP_CN3","PPP_CN4","PPP_CN5","PPP_EX1","PPP_EX2","PPP_EX3","PPP_EX4","PPP_EX5")]
dim(GenewiseCounts)


#trim_data_only_DPP
GenewiseCounts<- GenewiseCounts[c(1:13)]

head(GenewiseCounts)
#trim_outlier
#GenewiseCounts <-  GenewiseCounts[-c(2:4,7:9)] #outlier DPP_CN4, DPP_EX1, DPP_EX2,

head(GenewiseCounts)

##-------------------------------------------------------------------------
samplesheet <- read.delim("/Users/skailasam/Desktop/mnt/projects/Petunia_Tobacco_RNAseq_C3G-IM-427/01_Petunia/Analysis_report/samplesheet.tsv",
                          header = T,
                          sep = "\t",
                          row.names = 1
)
samplesheet$Sample_Name<- row.names(samplesheet)
#samplesheet$Sample_Plate <- NULL ;samplesheet$Sample_Well <- NULL;samplesheet$Description <- NULL; samplesheet$I7_Index_ID<- NULL; samplesheet$index<- NULL
#samplesheet$RepID<- c(rep(c(1:5), 4))
#colnames(samplesheet)<- c("SampleNameOld","SampleNameNew","RepID")
#samplesheet$NameID <- paste0(samplesheet$Sample_Name, samplesheet$RepID)
samplesheet
#trim_data_only_DPP
#samplesheet<- samplesheet[c(1:12),] 
#samplesheet<- samplesheet[-c(1:10,11:13,16,18,19),] # #outlier DPP_CN4
samplesheet
# out_file <- paste0(figdir,"/DataTables_summaryInfo.csv")
# sep <- ","
# cat("Note : Contains list of samples considered for this analysis after removing outliers\n", file = out_file)
# write.table(samplesheet, out_file, append = TRUE, sep = sep)

##-----Group---------------------------------------------------------------
anno <- data.frame(
  "Condition" = c("P1_1","P1_2","P1_3","P2_1","P2_2","P2_3","P2.5_1","P2.5_2","P2.5_3","P4_1","P4_2","P4_3"),
  #  "Condition" = c("DPP_CN","DPP_CN","DPP_CN","DPP_CN","DPP_EX","DPP_EX","DPP_EX","DPP_EX"), #outlier
  # "Condition" = c("DPP_CN","DPP_CN","DPP_EX","DPP_EX"), #outlier
  #  "Condition" = samplesheet$SampleNameNew,
  "Group" = samplesheet$Group,
  "RepID"     = samplesheet$RepID   #outlier
)
#trim samples
#anno <- anno[-c(2),]


group <- factor(anno$Group)
table(group)



## ----DGEList, message=FALSE----------------------------------------------
#Symbol  X06  X13   X4  X11   X5  X12   X10  X09  X08  X07  X03  X02  X01
#1       2     3    4    5    6   7     8     9    10   11  12    13  14

dge <-
  DGEList(GenewiseCounts[-c(1)], group = anno$Group, genes =
            GenewiseCounts[, 1, drop = FALSE])
#options(digits=3)
dge$samples
#rownames(dge$samples) <- c("S1","S2","S3","S4","S5","S6","S7","S8")
#delete duplicate genesymbols
o <- order(rowSums(dge$counts))
dge <- dge[o,]
d <- duplicated(dge$genes$Symbol)
dge <- dge[!d,]



## ----dropNAsymbols-------------------------------------------------------
dim(dge)
dge <- dge[!is.na(dge$genes), ]
dim(dge)

#trim samples
#dge <- DGEList(GenewiseCounts[c(2, 6, 7, 5, 3)], group=anno$group,genes=GenewiseCounts[,1,drop=FALSE])




## ----design---contrast----------------------------------------------------------

design1 = model.matrix( ~ 0 + Group, data = anno)
is.fullrank(design1)
design1

colnames(design1) <- make.names(colnames(design1))
contr.matrix <-
  makeContrasts(
    PhyII_vsI  = (GroupII  - GroupI),
    PhyIII_vsI  = (GroupIII - GroupI),
    PhyIV_vsI   = (GroupIV  - GroupI),
    PhyIII_vsII = (GroupIII - GroupII),
    PhyIV_vsII  = (GroupIV  - GroupII),
    PhyIV_vsIII = (GroupIV  - GroupIII),
    levels = colnames(design1)
  )
contr.matrix

# 
color4_custom<- c("#a6611a",
                  "#dfc27d",
                  "#80cdc1",
                  "#018571")
# The names argument tells the barplot to use the sample names on the x-axis
# The las argument rotates the axis names
cnt <- cnt + 1
png(file.path(
  figdir,
  paste0(figcount[cnt], "barplot_libray_size_countfilt1.png")
),width = 10, height = 8,units = 'in', res = 300)
col.cell <-
  c("#a6611a",
    "#dfc27d",
    "#80cdc1",
    "#018571")[anno$Group]
barplot(
  dge$samples$lib.size,
  names = anno$Condition,
  las = 2,
  col = col.cell
) +
  # Add a title to the plot
  title("Barplot of library sizes")
dev.off()




v <- voom(dge, design1  , plot = FALSE)

# PCA
##################################################################
exp = v$E
x = prcomp(as.matrix(t(exp)))
perc_pc <- x$sdev ^ 2 / sum(x$sdev ^ 2) * 100

pca.ap = x$x
pca.ap = data.frame(pca.ap)
pca.ap[, "Group"] <- anno$Group
pca.ap[, "RepID"] <- anno$RepID

i = 1
j = 2
cnt <- cnt + 1
png(file = file.path(
  figdir,
  paste0(figcount[cnt], "PCA_sample_group.png")
),width = 10, height = 8,units = 'in', res = 300)

ggplot(pca.ap, aes(
  x = eval(parse(text = paste0("PC", i))),
  y = eval(parse(text = paste0("PC", j))),
  colour = Group,
  #   shape = RepID
)) +
  geom_point(stroke = 1.5,
             size = 3,
             alpha = 1) +
  scale_color_manual(values = c("#a6611a",
                                "#dfc27d",
                                "#80cdc1",
                                "#018571"
  )) +
  
  # scale_color_manual(values = c("dodgerblue3","seagreen2","firebrick3","#E69F00", "darkmagenta", "#ea10c9")) +
  scale_shape_manual(values = c(1,2,3,4)) +
  geom_text(
    aes(label = RepID),
    size = 4,
    check_overlap = FALSE,
    nudge_y = 3
  ) +
  xlab(paste0("PC", i, ": ", round(perc_pc[[i]], 2), "%")) +
  ylab(paste0("PC", j, ": ", round(perc_pc[[j]], 2), "%")) +
  #theme_bw(base_size = 15)  +
  theme(axis.text.x = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 20, angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
  # Change the appearance and the orientation angle
  # of axis tick labels
  theme(axis.text.x = element_text(face="plain", color="black",
                                   size=14, angle=0),
        axis.text.y = element_text(face="plain", color="black",
                                   size=14, angle=90)) +
  
  
  
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.line.x= element_line(colour = "black"),
    axis.line.y.left = element_line(colour = "black")
  ) + theme(aspect.ratio = 1)
#+ theme(legend.position = "none")
#hc <- hclust(dist(as.matrix((t(exp)))), "ward.D2")

dev.off()




### MDS plot

mds <- plotMDS(dge)
toplot <-
  data.frame(Dim1 = mds$x,
             Dim2 = mds$y,
             RepID = anno$RepID,
             Group = anno$Group)
cnt <- cnt + 1
png(file = file.path(
  figdir,
  paste0(figcount[cnt], "MDS_sample_group.png")
),width = 10, height = 8,units = 'in', res = 300)

ggplot(toplot, aes(Dim1, Dim2, colour = Group)) +
  geom_text(
    aes(label = RepID),
    size = 4,
    check_overlap = FALSE,
    nudge_y = 0.05
  ) +
  geom_point(stroke = 1.5,
             size = 3,
             alpha = 1) +
  #  title("MDS plot") +
  xlab("Leading logFC dim 1") +
  ylab("Leading logFC dim 2") +
  scale_color_manual(values = c("#a6611a",
                                "#dfc27d",
                                "#80cdc1",
                                "#018571"
  ))     +
  #theme_bw(base_size = 15)  +
  theme(axis.text.x = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 20, angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
  # Change the appearance and the orientation angle
  # of axis tick labels
  theme(axis.text.x = element_text(face="plain", color="black",
                                   size=14, angle=0),
        axis.text.y = element_text(face="plain", color="black",
                                   size=14, angle=90)) +
  
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.line.x= element_line(colour = "black"),
    axis.line.y.left = element_line(colour = "black")
  ) + theme(aspect.ratio = 1)
#+ theme(legend.position = "none")
#hc <- hclust(dist(as.matrix((t(exp)))), "ward.D2")
#  theme(aspect.ratio = 1) + theme(legend.position = "none")
dev.off()

dim(dge)
#options(digits=3)
dge$samples

#delete duplicate genesymbols
o <- order(rowSums(dge$counts))
dge <- dge[o,]
d <- duplicated(dge$genes$Symbol)
dge <- dge[!d,]
dim(dge)


## ----dropNAsymbols-------------------------------------------------------
dge <- dge[!is.na(dge$genes), ]
dim(dge)

rownames(dge) <- dge$genes$Symbol
## ----keep----------------------------------------------------------------

#option 1
#keep <- rowSums(cpm(dge) > 0.5) >= 2
countsPerMillion <- cpm(dge)
summary(countsPerMillion)

countCheck <- countsPerMillion > 1
head(countCheck)
keep <- which(rowSums(countCheck) >= 2)
dge <- dge[keep,]
summary(cpm(dge)) #compare this to the original summary
#xs4321  `1654`
#option 2
keep <-
  filterByExpr(dge) #keeps rows that have worthwhile counts in a minumum number

table(keep)
## ----filter--------------------------------------------------------------
dim(dge)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dim(dge)


## ----norm----------------------------------------------------------------
dge <- calcNormFactors(dge,method="TMM")
dge$samples

# # Obtain CPMs
myCPM <- cpm(dge$counts, log = T)
# # Have a look at the output
head(myCPM)
write.csv(myCPM,
          file = file.path(tabdir, "Table_logCPM_matrix.csv"),
          quote = TRUE, row.names = T
)

cnt <- cnt + 1
png(file.path(
  figdir,
  paste0(figcount[cnt], "barplot_libray_size_countfilt2.png")
),width = 10, height = 8,units = 'in', res = 300)
col.cell <-
  c("#a6611a",
    "#dfc27d",
    "#80cdc1",
    "#018571")[anno$Group]
barplot(
  dge$samples$lib.size,
  names = anno$Condition,
  las = 2,
  col = col.cell
) +
  # Add a title to the plot
  title("Barplot of library sizes")
dev.off()



v <- voom(dge, design1  , plot = FALSE)

# PCA
##################################################################
exp = v$E
x = prcomp(as.matrix(t(exp)))
perc_pc <- x$sdev ^ 2 / sum(x$sdev ^ 2) * 100

pca.ap = x$x
pca.ap = data.frame(pca.ap)
pca.ap[, "Group"] <- anno$Group
pca.ap[, "RepID"] <- anno$RepID

i = 1
j = 2
cnt <- cnt + 1
png(file = file.path(
  figdir,
  paste0(figcount[cnt], "PCA_sample_group.png")
),width = 10, height = 8,units = 'in', res = 300)

ggplot(pca.ap, aes(
  x = eval(parse(text = paste0("PC", i))),
  y = eval(parse(text = paste0("PC", j))),
  colour = Group,
  #   shape = RepID
)) +
  geom_point(stroke = 1.5,
             size = 3,
             alpha = 1) +
  scale_color_manual(values = c("#a6611a",
                                "#dfc27d",
                                "#80cdc1",
                                "#018571"
  )) +
  
  # scale_color_manual(values = c("dodgerblue3","seagreen2","firebrick3","#E69F00", "darkmagenta", "#ea10c9")) +
  scale_shape_manual(values = c(1,2,3,4)) +
  geom_text(
    aes(label = RepID),
    size = 4,
    check_overlap = FALSE,
    nudge_y = 3
  ) +
  xlab(paste0("PC", i, ": ", round(perc_pc[[i]], 2), "%")) +
  ylab(paste0("PC", j, ": ", round(perc_pc[[j]], 2), "%")) +
  #theme_bw(base_size = 15)  +
  theme(axis.text.x = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 20, angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
  # Change the appearance and the orientation angle
  # of axis tick labels
  theme(axis.text.x = element_text(face="plain", color="black", 
                                   size=14, angle=0),
        axis.text.y = element_text(face="plain", color="black", 
                                   size=14, angle=90)) +
  
  
  
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.line.x= element_line(colour = "black"),
    axis.line.y.left = element_line(colour = "black")
  ) + theme(aspect.ratio = 1)
#+ theme(legend.position = "none")
#hc <- hclust(dist(as.matrix((t(exp)))), "ward.D2")

dev.off()




### MDS plot


mds <- plotMDS(dge)
toplot <-
  data.frame(Dim1 = mds$x,
             Dim2 = mds$y,
             RepID = anno$RepID,
             Group = anno$Group)
cnt <- cnt + 1
png(file = file.path(
  figdir,
  paste0(figcount[cnt], "MDS_sample_group.png")
),width = 10, height = 8,units = 'in', res = 300)

ggplot(toplot, aes(Dim1, Dim2, colour = Group)) +
  geom_text(
    aes(label = RepID),
    size = 4,
    check_overlap = FALSE,
    nudge_y = 0.05
  ) +
  geom_point(stroke = 1.5,
             size = 3,
             alpha = 1) +
  #  title("MDS plot") +
  xlab("Leading logFC dim 1") +
  ylab("Leading logFC dim 2") +
  scale_color_manual(values = c("#a6611a",
                                "#dfc27d",
                                "#80cdc1",
                                "#018571"
  ))     +
  #theme_bw(base_size = 15)  +
  theme(axis.text.x = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 20, angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
  # Change the appearance and the orientation angle
  # of axis tick labels
  theme(axis.text.x = element_text(face="plain", color="black", 
                                   size=14, angle=0),
        axis.text.y = element_text(face="plain", color="black", 
                                   size=14, angle=90)) +
  
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.line.x= element_line(colour = "black"),
    axis.line.y.left = element_line(colour = "black")
  ) + theme(aspect.ratio = 1)
#+ theme(legend.position = "none")
#hc <- hclust(dist(as.matrix((t(exp)))), "ward.D2")
#  theme(aspect.ratio = 1) + theme(legend.position = "none")
dev.off()





# DEG
dge <- estimateDisp(dge, design1)
cnt <- cnt + 1
png(file = file.path(
  figdir,
  paste0(figcount[cnt], "Dispersion_estimate.png")
),width = 10, height = 8,units = 'in', res = 300)
plotBCV(dge, cex = 0.5,main="Dispersion Estimate")
dev.off()

v <- voom(dge, design1, plot = FALSE)
fit <-lmFit(v,design1)
plotSA(fit, main="voom: Mean-variance trend")

vfit<- contrasts.fit(fit,contrasts = contr.matrix)
efit <- eBayes(vfit)
cnt <- cnt + 1
png(file = file.path(
  figdir,
  paste0(figcount[cnt], "Mean-variance_trend.png")
),width = 10, height = 8,units = 'in', res = 300)
plotSA(efit, main="Final model: Mean-variance trend")
dev.off()




# fit <-glmFit(dge,design1)
# 
# lrt <- glmLRT(fit,coef = 1)
# topTags(lrt)
# print(summary(decideTests(lrt, p.value = aPVALcutoff, lfc = log2FCcutval)))
# lrt <- glmLRT(fit,coef = 2)
# topTags(lrt)
# print(summary(decideTests(lrt, p.value = aPVALcutoff, lfc = log2FCcutval)))
# lrt <- glmLRT(fit,coef = 3)
# topTags(lrt)
# print(summary(decideTests(lrt, p.value = aPVALcutoff, lfc = log2FCcutval)))

summary.table<- as.data.frame.matrix(summary(decideTests(efit, p.value = aPVALcutoff, lfc = log2FCcutval)))
write.table(summary.table,file.path(figdir,"Table_summary_DEG_stat.tsv"))

newdf<- rbind(summary.table, data.table(PhyII_vsI=paste0("*FC >= ",FCcutoff,"*FDR >= ",aPVALcutoff)),fill=T)
#rownames(newdf)<- NULL


dfrownames <- as.data.frame(c(row.names(summary.table),""))
colnames(dfrownames)<- c("Gene Type")
df.combo <- cbind(dfrownames,newdf)

cnt <- cnt + 1
png(file = file.path(
  figdir,
  paste0(figcount[cnt], "Table_summary_DEG_stat.png")
),width = 10, height = 8,units = 'in', res = 300)
grid.table(df.combo)
dev.off()

###Res1####
#create contrast list
contrast_list <- c("PhyII_vsI", "PhyIII_vsI","PhyIV_vsI","PhyIII_vsII","PhyIV_vsII","PhyIV_vsIII")
my_contrast_list <- list()
#DEG RAW for  all 6 comparison

for(y in 1:6){ #Loop through the  six comparisons
  
  # RES1 
  res1<-
    data.frame(topTable(
      efit,
      coef = contrast_list[y],
      adjust = "BH",
      n = nrow(dge$counts)
      # sort.by = "none"
    ))[, c("Symbol", "logFC", "P.Value", "adj.P.Val")]
  
  
  colnames(res1) <- paste( colnames(res1),contrast_list[y], sep = ".")
  
  write.csv(
    res1,
    file = file.path(tabdir, paste("Table",contrast_list[y],"_DEG_TotalValidGenes_in_Genome.csv",sep="_")),
    quote = TRUE, row.names = F
  )
  
  
  write.table(
    rownames(res1),
    file = file.path(tabdir, "genenames_DEG_TotalValidGenes_in_Genome_GO_bg.txt"),
    col.names = F,
    row.names = F,
    quote = FALSE,
    sep = "\t"
  )
  my_contrast_list[[y]] <- res1
}
###Res2####
# RES2

for(y in 1:6){ #Loop through the 
  
  res2<-
    data.frame(topTable(
      efit,
      coef = contrast_list[y],
      adjust = "BH",
      n = nrow(dge$counts) , p.value = aPVALcutoff, lfc = log2FCcutval
      # sort.by = "none"
    ))[, c("Symbol", "logFC", "P.Value", "adj.P.Val")]
  
  colnames(res2) <- paste( colnames(res2), contrast_list[y],sep = ".")
  
  # TOTAL DEG 
  write.csv(
    res2,
    file = file.path(tabdir, paste("Table_deg_",contrast_list[y],".csv",sep = "")),
    quote = TRUE, row.names = F
  )
  
  write.table(
    rownames(res2),
    file = file.path(tabdir, paste("genenames4_GO_deg_",contrast_list[y],".txt",sep = "")),
    col.names = F,
    row.names = F,
    quote = FALSE,
    sep = "\t"
  )
  
  # UP REG
  
  write.table(
    res2[which(res2[2]> log2FCcutval &
                 res2[4]< aPVALcutoff), ],
    file = file.path(tabdir, paste("Table_upreg_",contrast_list[y],".txt",sep = "")),
    col.names = F,
    row.names = F,
    quote = FALSE,
    sep = "\t"
  )
  
  
  write.table(
    rownames(res2[which(res2[2]> log2FCcutval &
                          res2[4]< aPVALcutoff), ]),
    file = file.path(tabdir, paste("genename4_GO_upreg_",contrast_list[y],".txt",sep = "")),
    col.names = F,
    row.names = F,
    quote = FALSE,
    sep = "\t"
  )
  
  #DOWN REG
  
  write.table(
    res2[which(res2[2]< (-log2FCcutval) &
                 res2[4]< aPVALcutoff), ],
    file = file.path(tabdir, paste("Table_downreg_",contrast_list[y],".txt",sep = "")),
    col.names = F,
    row.names = F,
    quote = FALSE,
    sep = "\t"
  )
  write.table(
    rownames(res2[which(res2[2]< (-log2FCcutval) &
                          res2[4]< aPVALcutoff), ]),
    file = file.path(tabdir, paste("genename4_GO_downreg_",contrast_list[y],".txt",sep = "")),
    col.names = F,
    row.names = F,
    quote = FALSE,
    sep = "\t"
  )
  
}

##### heatmap####
##################################################################
exp <- v$E
colnames(exp) <- paste0(colnames(exp), "_", anno$Condition)


cnt <- cnt + 1
for (p in 1:6) {
  deg.genes <-
    read.table(file.path(tabdir, paste0("genenames4_GO_deg_", contrast_list[p], ".txt")), header = F)
  print(nrow(deg.genes))
  deg <-
    exp[rownames(exp) %in% deg.genes[, 1], grep("_", colnames(exp))]
  
  deg <- t(scale(t(as.matrix(deg))))
  mtscaled <- deg
  dim(mtscaled)
  hc = hclust((dist(mtscaled)), method = "ward.D")
  
  col = colorRampPalette(c("blue", "white", "red"))
  png(file.path(figdir, paste0(figcount[cnt], "heatmap_avgExpr_", contrast_list[p], ".png")),width = 10, height = 8,units = 'in', res = 300)
  hm <-
    heatmap.2(
      mtscaled,
      labCol = anno$IDs ,
      symm = F,
      scale = 'none',
      main = paste0("Heatmap: ", contrast_list[p]),
      ylab = paste0("FC > ", FCcutoff, " ;  adj.Pval <  ", aPVALcutoff),
      cexRow = 0.75,
      cexCol = 0.75,
      trace = "none",
      dendrogram = 'column',
      Rowv = as.dendrogram(hc),
      #Colv = T,
      col = col(24),
      symbreaks = T,
      srtCol = 45,
      labRow = FALSE
    )
  dev.off()
}

#####logCPMheatmap####

lcpm <- cpm(dge, log = T)
cnt <- cnt + 1
newcnt2 <- cnt

for (i in c(1:6)) {
  Cntrl_vs_Treatment <- topTreat(efit, coef = i, n = Inf, sort.by = "P")
  Cntrl_vs_Treatment.topgenes <- rownames(Cntrl_vs_Treatment)[1:25]
  var1 <- which(row.names(v$genes) %in% Cntrl_vs_Treatment.topgenes)
  #  mycol <- colorpanel(1000,"blue","white","red")
  col = colorRampPalette(c("blue", "white", "red"))
  deg <- lcpm[var1, ]
  deg <- t(scale(t(as.matrix(deg))))
  mtscaled <- deg
  dim(mtscaled)
  hc = hclust((dist(mtscaled)), method = "ward.D")
  
  #rownames(mtscaled) <- mapIds(org.Mm.eg.db, rownames(dge), column="ENTREZID", keytype="ENSEMBL") 
  png(file.path(
    figdir,
    paste0(figcount[newcnt2], "heatmap_top25_logCPM", colnames(efit$contrasts)[i], ".png")
  ),width = 10, height = 10,units = 'in', res = 300)
  heatmap.2(
    mtscaled,
    labCol = anno$IDs ,
    symm = F,
    main = paste0("Heatmap: ", colnames(efit$contrasts)[i]," FC > ", FCcutoff, " ;  adj.Pval <  ", aPVALcutoff),
    xlab = "logCPM HeatMap",
    #  title  = ,
    scale = 'none',
    cexRow = 0.75,
    cexCol = 0.75,
    trace = "none",
    dendrogram = 'none',
    Rowv = as.dendrogram(hc),
    Colv = F,
    col = col(24),
    symbreaks = T,
    srtCol = 45,
    margins = c(8, 10)
  )
  dev.off()
}



#, labRow = rownames(deg[sample(),]


#####volcano plot####
volcano_panel <- list()
for ( z in 1:6){
  # fig1
  # z =2
  gene_list = my_contrast_list[[z]]
  #gene_list[is.na(gene_list)] <- 1
  gene_list$color_flag <-
    as.factor(
      ifelse(
        gene_list[2]< -log2FCcutval  &
          gene_list[4]< aPVALcutoff,
        "Downreg",
        ifelse(
          gene_list[2]> log2FCcutval  &
            gene_list[4]< aPVALcutoff,
          "Upreg" ,
          "Non-DE"
        )
      )
    )
  
  colnames(gene_list)<- c("Symbol","logFC","P.Value","adj.P.Val","color_flag" )
  gene_list$gene_names <- gene_list$Symbol
  #gene_list[!gene_list$Gene_Symbol %in% rownames(deg[1:30,]), "gene_names"] <- ""
  
  #convertlofFC to FC
  
  
  volcano_panel[[z]] <-  ggplot(data = gene_list,
                                aes(
                                  x = logFC,
                                  y = -log10(adj.P.Val),
                                  colour = color_flag
                                )) +
    #, label = gene_names
    geom_point(size = 0.75) +
    # geom_text(size=2, color="black", hjust=0, vjust=3, aes(logFC,-log10(P.Value),label=gene_names)) +
    coord_cartesian(xlim = c(-10, 10), ylim = c(0, 9)) +
    labs(title=paste(contrast_list[z]), x ="Fold-change (log2)", y ="P-value (-log10)") +
    xlab("Fold-change (log2)") +
    ylab("P-value (-log10)") +
    scale_colour_manual(values = c("dodgerblue4", "gray", "firebrick4")) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black")
    ) +
    theme(aspect.ratio = 1) + theme(legend.title = element_blank())
  #  geom_text_repel(size = 2.5, color = "black")
  
}

##-----------------------------------------------------
m1 <- grid.arrange(grobs = volcano_panel, ncol = 2) ## display plot
ggsave(m1, file = file.path( figdir,paste0(figcount[cnt+1],"volcanoplot_combined.png")))

################
# load libraries
library(plyr)
library(reshape2)
library(UpSetR)

# make a list of input files to be read
filelist = list.files(pattern = "genenames4_GO_deg_",path = c("/Users/skailasam/Desktop/mnt/projects/Petunia_Tobacco_RNAseq_C3G-IM-427/01_Petunia/Analysis_report/output_FC_2.0_FDR_0.05/Tables"))

# make a 3 column table of listname,gene,1
res<-lapply(filelist, function(x){
  data.frame(
    set=x,
    geneID=as.character(read.table(x)[,1]),
    val=1)
})

res<-ldply(res)

# turn the 3 column long table into a wide
res1<-acast(res,geneID~set,value.var="val",fill=0) 

# force as dataframe
res1<-as.data.frame(res1)

# 1st column must be a name
res1$name=rownames(res1)

# rearrange columns
res2=res1[,c(ncol(res1),1:(ncol(res1)-1))]

pdf("gene_intersections.pdf")
#make the plot with 6 sets
upset(res2,nsets = 6)
dev.off()


#####Venn####


dt <- decideTests(efit, p.value = aPVALcutoff, lfc = log2FCcutval)
datout <- summary(dt)
write.fit(efit,dt,
          file = file.path(tabdir, "Table_summary_listofgenes_results.tsv"))
de.common <-
  # dt[which(dt[, c(1)] != 0 | dt[, c(2)] != 0 | dt[, c(3)] != 0), ]
  dt[which(dt[, c(1)] != 0) , ]
nrow(de.common)

pdf("output_FC2.0_Pval0.05/loose.pdf")

vennDiagram(
  # de.common[, 1:3], # edit based on no. of sample
  de.common[, 1:5],
  include = c("up", "down"),
  cex = 1,
  lwd = 2,
  counts.col = c("red", "blue"),
  circle.col = c("#f4a582", "#99d594", "#92c5de", "#0571b0")
)
dev.off()



#panel.ul <- image_read_pdf("strict.pdf")
panel.ur <- image_read_pdf("output_FC2.0_Pval0.05/loose.pdf")

whole <- c(panel.ur) %>%
  image_append()

#pdf(file.path(figdir,paste0(figcount[cnt],"venndiagram.pdf"))) # saves the current plotting device to a pdf file.
#plot(whole) # plot the external object generated in ImageMagick to R's plotting device
#dev.off()
#file.remove("output_FC2.0_Pval0.05/loose.pdf")

combo2 <- ggdraw() +
  draw_label(
    "Figure: venndiagram  ",
    fontface = 'bold',
    x = 0.1,
    y = 1,
    hjust = 0,
    vjust = 1
  ) +
  draw_label(
    "(A )",
    x = 0.1,
    y = .9,
    hjust = 0,
    vjust = 1
  ) +
  draw_label(
    paste0("FC > ", FCcutoff, " ;  adj.Pval <  ", aPVALcutoff),
    x = 0.1,
    y = .85,
    hjust = 0,
    vjust = 1
  ) +
  draw_image(
    whole,
    x = 0,
    y = 0,
    width = 1,
    height = 0.75
  )
cnt <- cnt + 1
save_plot(file.path(figdir, paste0(figcount[cnt], "venndiagram_combined.png")),
          combo2,
          base_width = 10,
          base_height = 8)




