Differentially Expressed genes using egdeR  Document:Use the egdeR user guide

#Provided we have a count.txt and design.txt 

Packages 
install biocMananger
install edgeR and load library(edgeR)
install pheatmap and load

Commands

 counts <- read.table("counts.txt",header="TRUE",sep="\t",row.names=1)
 design <-read.table("design.txt",header="TRUE",sep="\t",row.names=1)
 group <-factor(design$Group)
 y <- DGEList(counts=counts,group=group)
 keep <- filterByExpr(y)
 y <- y[keep,keep.lib.sizes=FALSE]
 y <- normLibSizes(y)
 plotMDS(y)

 design <- model.matrix(~group)
 y <- estimateDisp(y, design)
 plotBCV(y)
 fit <- glmQLFit(y,design)
 qlf <- glmQLFTest(fit,coef=2)
 topTags(qlf)
 summary(decideTests(qlf))
 plotMD(qlf) 

 top_genes <- topTags(qlf, n = 30)
 top_gene_names <- rownames(top_genes$table)
 cpm_data <- cpm(y, log = TRUE,prior.count=1)
 cpm_top <- cpm_data[top_gene_names, ]
 color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
> pheatmap(cpm_top, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         main = "Heatmap of Top 30 DEGs", 
         scale = "row",
         color = color_palette,
	  )
pheatmap(cpm_top,
 cluster_rows= TRUE,
 cluster_cols= TRUE,
 show_rownames= TRUE,
 show_colnames= TRUE,
 color=color_palette,
 main="deg of 30 genes",
 annotation_col= annotation_col,
 clustering_distance_cols="euclidean",
 clustering_method="complete")