Treeplot and Treemap

ego <- enrichGO(gene          = gene.df$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
        readable      = TRUE)


 ego<-pairwise_termsim(ego)
 ego_sim_df<-as.data.frame(ego@result)

Treemap

library(treemap)
 treemap(ego_sim_df,
+ index="Description",
+ vSize="Count",
+ vColor="p.adjust",
+ title="GO Term Enrichemnt",
+ palette="Blues")

Treeplot

treeplot(ego)

