Gene Ontology Enrichment Analysis 

#Provided a gene list with the logFC score
 
Packages
install and load clusterProfiler and org.Hs.eg.db from BiocManager
install and load readr

Commands
 
 df<-read.csv("geos.csv")
 head(df)
gene.df <- bitr(df$gene_symbol, fromType = "SYMBOL",   
    toType = ("ENTREZID"),
    OrgDb = org.Hs.eg.db)

Gene Ontology Classification
 ggo <- groupGO(gene     = gene.df$ENTREZID,
+                OrgDb    = org.Hs.eg.db,
+                ont      = "CC",
+                level    = 3,
+                readable = TRUE)
> head(ggo)

Gene Over representation Analysis
 ego2 <- enrichGO(gene         = gene.df$ENTREZID,
+                 OrgDb         = org.Hs.eg.db,
+                 keyType       = 'ENTREZID',
+                 ont           = "CC",
+                 pvalueCutoff  = 0.01,
+                 pAdjustMethod = "BH",
+                 qvalueCutoff  = 0.05)
> head(ego2, 3)
dotplot(ego2)

Gene Enrichment Analysis

create a genelist
id <- bitr(gene, fromType = "ENTREZID", 
      toType = ("ENTREZID"),
      OrgDb = org.Hs.eg.db)
geneList<-setNames(df_merged$logFC,df_merged$ENTREZID)
sorted_geneList<-sort(geneList,decreasing=TRUE)
ego3 <- gseGO(geneList     = sorted_geneList,
              OrgDb        = org.Hs.eg.db,
              ont          = "CC",
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)
dotplot(ego3)
