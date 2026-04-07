Heatmap

## load edgeR
library('edgeR')

### Input raw count data into edgeR ###
rawdata <- read.table("counts.txt", header=TRUE, row.names=1)

dim(rawdata)

## Read a design file with sample information ###
design<- read.table("design.txt", sep ="\t", row.names=1, header =TRUE)

### Make class label ###
group <- factor(design$Group)

# Make DGEList object ###
y <- DGEList(counts=rawdata, group=group)

# Transformation form raw scale ###
### Before and after pre-processing ###
cpm <- cpm(y)
lcpm <- cpm(y, log=TRUE)

### Filter out lowly expressed genes ###
keep <- rowSums(cpm(y)>2) >= 3
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]

########################################

### Normalization for composition bias
y <- calcNormFactors(y)
y$samples

##### Plot MDS
plotMDS(y, col = as.numeric(y$samples$group), cex = 1, labels=NULL)

############### Estimate dispersion ###
###y <- estimateDisp(y,deisgn)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)

### Result for control - test ####
et1 <- exactTest(y, pair=c("control","test"))
res <- topTags(et1)

con.test <- topTags(et1, n=nrow(et1$table))
write.table(con.test, file="controlVStest.txt", sep = "\t")

#### Heatmap clustering ####
logCPM <- cpm(y, prior.count=2, log=TRUE)
#rownames(logCPM) <- y$genes$Symbol
#colnames(logCPM) <- paste(y$samples$group, 1:2, sep="-")

o1 <- order(et1$table$PValue)
logCPM.con.test <- logCPM[o1[1:30],]

logCPM.con.test <- t(scale(t(logCPM.con.test)))
library(gplots)
col.pan <- colorpanel(100, "yellow", "white", "pink")
pdf("con_D7COC_heatmap.pdf")
heatmap.2(logCPM.con.test, col=col.pan, Rowv=TRUE, scale="none", trace="none", dendrogram="both", cexRow=1, cexCol=1.4, density.info="none", margin=c(10,9), lhei=c(2,10), lwid=c(2,6))
dev.off()
