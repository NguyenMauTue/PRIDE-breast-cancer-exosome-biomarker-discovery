############################################################
# 13 Driver score via PCA
############################################################

pca_testing = data.frame(
  logFC = limmanetwork_df$logFC,
  nlog10FDR = -log10(limmanetwork_df$adj.P.Val),
  degree = log(limmanetwork_df$degree + 1),
  betweenness = log(limmanetwork_df$betweenness + 1),
  row.names = limmanetwork_df$Symbol
)

pca_testing_data = prcomp(
  pca_testing,
  scale = TRUE
)

pca_testing_x =
  as.data.frame(pca_testing_data$x)

var =
  pca_testing_data$sdev^2 /
  sum(pca_testing_data$sdev^2)


############################################################
# Driver score
############################################################

alpha = var[2] / var[1]

DriverScore =
  pca_testing_x$PC1 +
  (alpha * pca_testing_x$PC2)
