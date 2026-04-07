Pathway Enrichment Analysis 

#Provided a gene list with the logFC score
 
Packages

install and load clusterProfiler and org.Hs.eg.db from BiocManager
install and load readr

Commands

geneList<-setNames(df_merged$logFC,df_merged$ENTREZID)
sorted_geneList<-sort(geneList,decreasing=TRUE)

KEGG ORA

kk <- enrichKEGG(gene         = id$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
head(kk)
dotplot(kk)

KEGG Gene set enrichment

kk2 <- gseKEGG(geneList     = sorted_geneList,
               organism     = 'hsa',
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
head(kk2)
dotplot(kk2)

for viewing a pathway with pathview
hsa05200 <- pathview(gene.data  = geneList,
                     pathway.id = "hsa05200",
                     species    = "hsa",
                     limit      = list(gene=max(abs(geneList)), cpd=1))

Wikipathways ORA

 data(geneList, package="DOSE")
 gene <- names(geneList)[abs(geneList) > 2]
 enrichWP(gene, organism = "Homo sapiens") 


Wikipathways Geneset enrichment Analysis

 gseWP(geneList, organism = "Homo sapiens")
