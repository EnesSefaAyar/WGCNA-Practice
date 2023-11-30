# Identify Clusters Custom

```{r}
best_cluster <- function(orderedCorrMatrix, maxClusterSize = 100) {
  starts <- c()
  ends <- c()
  for (t in 1:5) {
    scores <- c()
    for (i in c(seq(from = 1, to = nrow(orderedCorrMatrix) - maxClusterSize, by = 1))) {
      for (j in c(seq(from = 1, to = maxClusterSize, by = 1))) {
        score <- mean(abs(orderedCorrMatrix[i:(i+j), i:(i+j)]))*(j^0.5)
        scores <- c(scores, score)
      }
    }
    max_i <- which.max(scores) %/% maxClusterSize
    max_j <- (length(scores) - max_i*(nrow(orderedCorrMatrix) - 
                                        maxClusterSize)) %% maxClusterSize 
    starts <- c(starts, max_i)
    ends <-  c(ends, max_i + max_j)
    
    orderedCorrMatrix <- orderedCorrMatrix[-max_i:-(max_i + max_j), -max_i:-(max_i + max_j)]
  }
  
  df <- data.frame(starts, ends)
  return(df)
}

indices <- best_cluster(ordered_cor_df)

indices[3, 2] - indices[3, 1]

cluster_members_CPM <- rownames(ordered_cor_df)[indices[3, 1]:indices[3, 2]]

# Map to Gene names
up <- UniProt.ws(taxId = 9606)
res <- select(
  x = up,
  keys = cluster_members_CPM,
  to = "Gene_Name"
)

# Subplot of module correlation
temp_df <- ordered_cor_df[cluster_members_CPM, cluster_members_CPM]
if (length(res$To) == length(rownames(temp_df))){
  rownames(temp_df) <- res$To 
  colnames(temp_df) <- res$To 
}
pheatmap(temp_df,
         color = colorRampPalette(c("#4d72af", "white", "#c44e53")) (101),
         breaks = seq(from = -1, to = 1, length.out = 101),
         border_color = NA,legend = FALSE, show_colnames = TRUE, 
         show_rownames = TRUE,
         cluster_rows = FALSE, cluster_cols = FALSE,
         cellwidth = 500/nrow(temp_df), cellheight = 500/nrow(temp_df),
         fontfamily = 'sans', fontsize = 10,
         display_numbers = round(temp_df, 2))
```
