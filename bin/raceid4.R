#RaceID4 
library(tidyverse)
library(viridis)
library(RaceID)
library(Matrix)
library(data.table)
library(Seurat)

date = Sys.Date()
#download counts file from url: https://stackoverflow.com/questions/28986150/downloading-and-extracting-gz-data-file-using-r
#aut_counts dataset

#url <- "https://cells.ucsc.edu/autism/rawMatrix.zip"

#counts <- Read10X("data/")
counts <- readMM("data/filtered_count_matrix.mtx")
rownames(counts) <- readLines("data/filtered_gene_row_names.txt")
metadata <- read.delim("data/filtered_column_metadata2.txt")

counts <- counts[, which(metadata$amyloid.group == "low")]

which(genes == "TMEM119")
hist(counts["TMEM119",])
hist(counts["SLC2A5",])
hist(counts["P2RY12",])
sum(counts["TMEM119",]>0)
sum(counts["SLC2A5",]>0)
sum(counts["P2RY12",]>0)


micr_counts <- counts[,which(counts["SLC2A5",]>0 | counts["P2RY12",]>0 | counts["TMEM119",]>0)]


save(micr_counts, file = "data/mathys_ad_nuc_seq-microglia.Robj")

load("data/mathys_ad_nuc_seq-microglia.Robj")

prdata <- as.data.frame(as.matrix(micr_counts))

sc <- SCseq(prdata)
# filtering of expression data
a <- apply(prdata, 2, sum)
sc <- filterdata(sc, mintotal=quantile(a, 0.1)) # exlcude the lower quartile of the cells
sc <- CCcorrect(sc, 
                dimR = T, 
                nComp = 20,
                CGenes = c('JUN',
                           'FOS',
                           'ZFP36',
                           'HSPA1A|HSPA1B',
                           'DUSP1',
                           'EGR1',
                           'MALAT1'))

sc <- compdist(sc,metric="pearson")
sc <- clustexp(sc) 

plotsaturation(sc,disp=FALSE)
plotsaturation(sc,disp=TRUE)
plotjaccard(sc)

sc <- clustexp(sc,cln=14,sat=FALSE) 
sc <- findoutliers(sc)
plotbackground(sc)
plotsensitivity(sc)
plotoutlierprobs(sc)
ord_clust <- clustheatmap(sc)
save(ord_clust, file = 'data/ord_clust.Robj')

pdf(paste0('plots/heatmaps/clustheatmap.pdf'))
clustheatmap(sc, final = T)
dev.off()

sc <- comptsne(sc)
sc <- compfr(sc,knn=10)

plotmap(sc)
plotmap(sc,fr=TRUE)
dev.off()

plotexpmap(sc,"MRC1",logsc=F,fr=F)
plotexpmap(sc,"LYVE1",logsc=F,fr=F)
plotexpmap(sc,"CD163",logsc=F,fr=F)
plotexpmap(sc,"TMEM119",logsc=F,fr=F)
plotexpmap(sc,"CX3CR1",logsc=F,fr=F)
plotexpmap(sc,"PTPRC",logsc=F,fr=F)
plotexpmap(sc,"CD3E",logsc=F,fr=F)
plotexpmap(sc,"ITGAM",logsc=F,fr=F)
plotexpmap(sc,"CD8A",logsc=F,fr=F)
plotexpmap(sc,"CD4",logsc=F,fr=F)
plotexpmap(sc,"P2RY12",logsc=F,fr=F)
plotexpmap(sc,"SLC2A5",logsc=F,fr=F)
plotexpmap(sc,"^EGR1",logsc=F,fr=F)
plotexpmap(sc,"JUN",logsc=F,fr=F)
plotexpmap(sc,"GPR34",logsc=F,fr=F)

dg <- clustdiffgenes(sc,4,pvalue=.01)
head(dg,25)
types <- sub("(\\_\\d+)$","", colnames(sc@ndata))
genes <- head(rownames(dg)[dg$fc>1],10)
plotmarkergenes(sc,genes,samples=types)

#Save sc file
save(sc, file = 'data/sc.Robj')

micr_ids <- names(sc@cpart)[sc@cpart %in% c(8:13)]
write_csv(as.data.frame(micr_ids), "data/microglia-cell-ids.csv")
