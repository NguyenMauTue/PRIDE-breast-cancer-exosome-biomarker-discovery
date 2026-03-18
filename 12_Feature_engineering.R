############################################################
# 12 Feature engineering
############################################################

limmanetwork_df =
  read.csv("../results/limma_network_table.csv")


exnetscore_df =
  limmanetwork_df[
    ,c("UNIPROT","logFC","degree","betweenness")
  ]


exnetscore_df$ExprScore =
  as.vector(
    abs(scale(exnetscore_df$logFC))
  )


exnetscore_df$NetScore =
  as.vector(
    scale(log(exnetscore_df$degree+1)) +
    scale(log(exnetscore_df$betweenness+1))
  )


write.csv(
  exnetscore_df,
  "../results/expression_network_features.csv",
  row.names=FALSE
)
