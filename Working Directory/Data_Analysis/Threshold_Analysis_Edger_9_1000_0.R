# Step 1: Call Library (edgeR)
library(edgeR)

# Step 2: Load Dataset
count_data <- read.table("RAW data/9_1000_0.tsv", header=TRUE, row.names=1)
groups <- factor(rep(1:2, each=9))

# Step 3: Create DGEList Object
y <- DGEList(counts=count_data, group=groups)

# Step 4: Filter out low count genes
keep <- filterByExpr(y)
y <- y[keep, ]

# Step 5: Normalize Values
y <- calcNormFactors(y)

# Step 6: Initialize Matrix for Values
design <- model.matrix(~groups)

# Step 7: EdgeR Dispersion test
y <- estimateDisp(y, design)

# Step 8: Extract Differential Expression data
et <- exactTest(y)

# Step 9: Summarize with DecideTestsDGE
p_value <- 0.05
decided <- decideTestsDGE(et, adjust.method="BH", p.value=p_value)
summary(decided)

# Step 10: Print top genes
topTags(et)

# Step 11: Load Metadata
meta_data <- read.table("RAW data/9_1000_0_meta.tsv", header=TRUE, row.names=1)
edgeR_results <- topTags(et, n=Inf)
annotated_results <- merge(edgeR_results, meta_data, by="row.names", all.x=TRUE)
rownames(annotated_results) <- annotated_results$Row.names
annotated_results$Row.names <- NULL

# Step 12: Comparison using LogFC
# Upregulated
detected_up_edgeR <- rownames(annotated_results[annotated_results$logFC > 0 & annotated_results$FDR < p_value,])
meta_up <- rownames(annotated_results[annotated_results$upregulation == 1,])
common_up <- intersect(detected_up_edgeR, meta_up)
# Downregulated
detected_down_edgeR <- rownames(annotated_results[annotated_results$logFC < 0 & annotated_results$FDR < p_value,])
meta_down <- rownames(annotated_results[annotated_results$downregulation == 1,])
common_down <- intersect(detected_down_edgeR, meta_down)

# Step 14: Accuracy and Precision Matrix
true_positives <- length(common_up) + length(common_down)
false_positives <- length(setdiff(detected_up_edgeR, meta_up)) + length(setdiff(detected_down_edgeR, meta_down))
true_negatives <- nrow(meta_data) - (length(meta_up) + length(meta_down)) - false_positives
false_negatives <- length(setdiff(meta_up, detected_up_edgeR)) + length(setdiff(meta_down, detected_down_edgeR))
accuracy <- (true_positives + true_negatives) / nrow(meta_data)
precision <- true_positives / (true_positives + false_positives)
recall <- true_positives / (true_positives + false_negatives)
f1_score <- 2 * ((precision * recall) / (precision + recall))
FDR <- false_positives/(true_positives+false_positives) ## Precision = 1- FDR 

# Step 15: Output Metrics to CSV for further analysis outside R
metrics_df_1 <- data.frame(
  Threshold = p_value,
  True_Positives = true_positives,
  False_Positives = false_positives,
  True_Negatives = true_negatives,
  False_Negatives = false_negatives,
  Accuracy = accuracy,
  Precision = precision,
  Recall = recall,
  F1_Score = f1_score,
  FDR = FDR,
  Experiment = "9_1000_0"
)

#p_value <- 0.04
p_value <- 0.04
decided <- decideTestsDGE(et, adjust.method="BH", p.value=p_value)
summary(decided)
edgeR_results <- topTags(et, n=Inf)
annotated_results <- merge(edgeR_results, meta_data, by="row.names", all.x=TRUE)
rownames(annotated_results) <- annotated_results$Row.names
annotated_results$Row.names <- NULL

# Upregulated
detected_up_edgeR <- rownames(annotated_results[annotated_results$logFC > 0 & annotated_results$FDR < p_value,])
meta_up <- rownames(annotated_results[annotated_results$upregulation == 1,])
common_up <- intersect(detected_up_edgeR, meta_up)
# Downregulated
detected_down_edgeR <- rownames(annotated_results[annotated_results$logFC < 0 & annotated_results$FDR < p_value,])
meta_down <- rownames(annotated_results[annotated_results$downregulation == 1,])
common_down <- intersect(detected_down_edgeR, meta_down)

true_positives <- length(common_up) + length(common_down)
false_positives <- length(setdiff(detected_up_edgeR, meta_up)) + length(setdiff(detected_down_edgeR, meta_down))
true_negatives <- nrow(meta_data) - (length(meta_up) + length(meta_down)) - false_positives
false_negatives <- length(setdiff(meta_up, detected_up_edgeR)) + length(setdiff(meta_down, detected_down_edgeR))
accuracy <- (true_positives + true_negatives) / nrow(meta_data)
precision <- true_positives / (true_positives + false_positives)
recall <- true_positives / (true_positives + false_negatives)
f1_score <- 2 * ((precision * recall) / (precision + recall))
FDR <- false_positives/(true_positives+false_positives) ## Precision = 1- FDR 

metrics_df_2 <- data.frame(
  Threshold = p_value,
  True_Positives = true_positives,
  False_Positives = false_positives,
  True_Negatives = true_negatives,
  False_Negatives = false_negatives,
  Accuracy = accuracy,
  Precision = precision,
  Recall = recall,
  F1_Score = f1_score,
  FDR = FDR,
  Experiment = "9_1000_0"
)

#p_value <- 0.03
p_value <- 0.03
decided <- decideTestsDGE(et, adjust.method="BH", p.value=p_value)
summary(decided)
edgeR_results <- topTags(et, n=Inf)
annotated_results <- merge(edgeR_results, meta_data, by="row.names", all.x=TRUE)
rownames(annotated_results) <- annotated_results$Row.names
annotated_results$Row.names <- NULL

# Upregulated
detected_up_edgeR <- rownames(annotated_results[annotated_results$logFC > 0 & annotated_results$FDR < p_value,])
meta_up <- rownames(annotated_results[annotated_results$upregulation == 1,])
common_up <- intersect(detected_up_edgeR, meta_up)
# Downregulated
detected_down_edgeR <- rownames(annotated_results[annotated_results$logFC < 0 & annotated_results$FDR < p_value,])
meta_down <- rownames(annotated_results[annotated_results$downregulation == 1,])
common_down <- intersect(detected_down_edgeR, meta_down)

true_positives <- length(common_up) + length(common_down)
false_positives <- length(setdiff(detected_up_edgeR, meta_up)) + length(setdiff(detected_down_edgeR, meta_down))
true_negatives <- nrow(meta_data) - (length(meta_up) + length(meta_down)) - false_positives
false_negatives <- length(setdiff(meta_up, detected_up_edgeR)) + length(setdiff(meta_down, detected_down_edgeR))
accuracy <- (true_positives + true_negatives) / nrow(meta_data)
precision <- true_positives / (true_positives + false_positives)
recall <- true_positives / (true_positives + false_negatives)
f1_score <- 2 * ((precision * recall) / (precision + recall))
FDR <- false_positives/(true_positives+false_positives) ## Precision = 1- FDR 

metrics_df_3 <- data.frame(
  Threshold = p_value,
  True_Positives = true_positives,
  False_Positives = false_positives,
  True_Negatives = true_negatives,
  False_Negatives = false_negatives,
  Accuracy = accuracy,
  Precision = precision,
  Recall = recall,
  F1_Score = f1_score,
  FDR = FDR,
  Experiment = "9_1000_0"
)

#p_value <- 0.02
p_value <- 0.02
decided <- decideTestsDGE(et, adjust.method="BH", p.value=p_value)
summary(decided)
edgeR_results <- topTags(et, n=Inf)
annotated_results <- merge(edgeR_results, meta_data, by="row.names", all.x=TRUE)
rownames(annotated_results) <- annotated_results$Row.names
annotated_results$Row.names <- NULL

# Upregulated
detected_up_edgeR <- rownames(annotated_results[annotated_results$logFC > 0 & annotated_results$FDR < p_value,])
meta_up <- rownames(annotated_results[annotated_results$upregulation == 1,])
common_up <- intersect(detected_up_edgeR, meta_up)
# Downregulated
detected_down_edgeR <- rownames(annotated_results[annotated_results$logFC < 0 & annotated_results$FDR < p_value,])
meta_down <- rownames(annotated_results[annotated_results$downregulation == 1,])
common_down <- intersect(detected_down_edgeR, meta_down)

true_positives <- length(common_up) + length(common_down)
false_positives <- length(setdiff(detected_up_edgeR, meta_up)) + length(setdiff(detected_down_edgeR, meta_down))
true_negatives <- nrow(meta_data) - (length(meta_up) + length(meta_down)) - false_positives
false_negatives <- length(setdiff(meta_up, detected_up_edgeR)) + length(setdiff(meta_down, detected_down_edgeR))
accuracy <- (true_positives + true_negatives) / nrow(meta_data)
precision <- true_positives / (true_positives + false_positives)
recall <- true_positives / (true_positives + false_negatives)
f1_score <- 2 * ((precision * recall) / (precision + recall))
FDR <- false_positives/(true_positives+false_positives) ## Precision = 1- FDR 

metrics_df_4 <- data.frame(
  Threshold = p_value,
  True_Positives = true_positives,
  False_Positives = false_positives,
  True_Negatives = true_negatives,
  False_Negatives = false_negatives,
  Accuracy = accuracy,
  Precision = precision,
  Recall = recall,
  F1_Score = f1_score,
  FDR = FDR,
  Experiment = "9_1000_0"
)

#p_value <- 0.06
p_value <- 0.06
decided <- decideTestsDGE(et, adjust.method="BH", p.value=p_value)
summary(decided)
edgeR_results <- topTags(et, n=Inf)
annotated_results <- merge(edgeR_results, meta_data, by="row.names", all.x=TRUE)
rownames(annotated_results) <- annotated_results$Row.names
annotated_results$Row.names <- NULL

# Upregulated
detected_up_edgeR <- rownames(annotated_results[annotated_results$logFC > 0 & annotated_results$FDR < p_value,])
meta_up <- rownames(annotated_results[annotated_results$upregulation == 1,])
common_up <- intersect(detected_up_edgeR, meta_up)
# Downregulated
detected_down_edgeR <- rownames(annotated_results[annotated_results$logFC < 0 & annotated_results$FDR < p_value,])
meta_down <- rownames(annotated_results[annotated_results$downregulation == 1,])
common_down <- intersect(detected_down_edgeR, meta_down)

true_positives <- length(common_up) + length(common_down)
false_positives <- length(setdiff(detected_up_edgeR, meta_up)) + length(setdiff(detected_down_edgeR, meta_down))
true_negatives <- nrow(meta_data) - (length(meta_up) + length(meta_down)) - false_positives
false_negatives <- length(setdiff(meta_up, detected_up_edgeR)) + length(setdiff(meta_down, detected_down_edgeR))
accuracy <- (true_positives + true_negatives) / nrow(meta_data)
precision <- true_positives / (true_positives + false_positives)
recall <- true_positives / (true_positives + false_negatives)
f1_score <- 2 * ((precision * recall) / (precision + recall))
FDR <- false_positives/(true_positives+false_positives) ## Precision = 1- FDR 

metrics_df_5 <- data.frame(
  Threshold = p_value,
  True_Positives = true_positives,
  False_Positives = false_positives,
  True_Negatives = true_negatives,
  False_Negatives = false_negatives,
  Accuracy = accuracy,
  Precision = precision,
  Recall = recall,
  F1_Score = f1_score,
  FDR = FDR,
  Experiment = "9_1000_0"
)

#p_value <- 0.07
p_value <- 0.07
decided <- decideTestsDGE(et, adjust.method="BH", p.value=p_value)
summary(decided)
edgeR_results <- topTags(et, n=Inf)
annotated_results <- merge(edgeR_results, meta_data, by="row.names", all.x=TRUE)
rownames(annotated_results) <- annotated_results$Row.names
annotated_results$Row.names <- NULL

# Upregulated
detected_up_edgeR <- rownames(annotated_results[annotated_results$logFC > 0 & annotated_results$FDR < p_value,])
meta_up <- rownames(annotated_results[annotated_results$upregulation == 1,])
common_up <- intersect(detected_up_edgeR, meta_up)
# Downregulated
detected_down_edgeR <- rownames(annotated_results[annotated_results$logFC < 0 & annotated_results$FDR < p_value,])
meta_down <- rownames(annotated_results[annotated_results$downregulation == 1,])
common_down <- intersect(detected_down_edgeR, meta_down)

true_positives <- length(common_up) + length(common_down)
false_positives <- length(setdiff(detected_up_edgeR, meta_up)) + length(setdiff(detected_down_edgeR, meta_down))
true_negatives <- nrow(meta_data) - (length(meta_up) + length(meta_down)) - false_positives
false_negatives <- length(setdiff(meta_up, detected_up_edgeR)) + length(setdiff(meta_down, detected_down_edgeR))
accuracy <- (true_positives + true_negatives) / nrow(meta_data)
precision <- true_positives / (true_positives + false_positives)
recall <- true_positives / (true_positives + false_negatives)
f1_score <- 2 * ((precision * recall) / (precision + recall))
FDR <- false_positives/(true_positives+false_positives) ## Precision = 1- FDR 

metrics_df_6 <- data.frame(
  Threshold = p_value,
  True_Positives = true_positives,
  False_Positives = false_positives,
  True_Negatives = true_negatives,
  False_Negatives = false_negatives,
  Accuracy = accuracy,
  Precision = precision,
  Recall = recall,
  F1_Score = f1_score,
  FDR = FDR,
  Experiment = "9_1000_0"
)

#p_value <- 0.08
p_value <- 0.08
decided <- decideTestsDGE(et, adjust.method="BH", p.value=p_value)
summary(decided)
edgeR_results <- topTags(et, n=Inf)
annotated_results <- merge(edgeR_results, meta_data, by="row.names", all.x=TRUE)
rownames(annotated_results) <- annotated_results$Row.names
annotated_results$Row.names <- NULL

# Upregulated
detected_up_edgeR <- rownames(annotated_results[annotated_results$logFC > 0 & annotated_results$FDR < p_value,])
meta_up <- rownames(annotated_results[annotated_results$upregulation == 1,])
common_up <- intersect(detected_up_edgeR, meta_up)
# Downregulated
detected_down_edgeR <- rownames(annotated_results[annotated_results$logFC < 0 & annotated_results$FDR < p_value,])
meta_down <- rownames(annotated_results[annotated_results$downregulation == 1,])
common_down <- intersect(detected_down_edgeR, meta_down)

true_positives <- length(common_up) + length(common_down)
false_positives <- length(setdiff(detected_up_edgeR, meta_up)) + length(setdiff(detected_down_edgeR, meta_down))
true_negatives <- nrow(meta_data) - (length(meta_up) + length(meta_down)) - false_positives
false_negatives <- length(setdiff(meta_up, detected_up_edgeR)) + length(setdiff(meta_down, detected_down_edgeR))
accuracy <- (true_positives + true_negatives) / nrow(meta_data)
precision <- true_positives / (true_positives + false_positives)
recall <- true_positives / (true_positives + false_negatives)
f1_score <- 2 * ((precision * recall) / (precision + recall))
FDR <- false_positives/(true_positives+false_positives) ## Precision = 1- FDR 

metrics_df_7 <- data.frame(
  Threshold = p_value,
  True_Positives = true_positives,
  False_Positives = false_positives,
  True_Negatives = true_negatives,
  False_Negatives = false_negatives,
  Accuracy = accuracy,
  Precision = precision,
  Recall = recall,
  F1_Score = f1_score,
  FDR = FDR,
  Experiment = "9_1000_0"
)

#p_value <- 0.09
p_value <- 0.09
decided <- decideTestsDGE(et, adjust.method="BH", p.value=p_value)
summary(decided)
edgeR_results <- topTags(et, n=Inf)
annotated_results <- merge(edgeR_results, meta_data, by="row.names", all.x=TRUE)
rownames(annotated_results) <- annotated_results$Row.names
annotated_results$Row.names <- NULL

# Upregulated
detected_up_edgeR <- rownames(annotated_results[annotated_results$logFC > 0 & annotated_results$FDR < p_value,])
meta_up <- rownames(annotated_results[annotated_results$upregulation == 1,])
common_up <- intersect(detected_up_edgeR, meta_up)
# Downregulated
detected_down_edgeR <- rownames(annotated_results[annotated_results$logFC < 0 & annotated_results$FDR < p_value,])
meta_down <- rownames(annotated_results[annotated_results$downregulation == 1,])
common_down <- intersect(detected_down_edgeR, meta_down)

true_positives <- length(common_up) + length(common_down)
false_positives <- length(setdiff(detected_up_edgeR, meta_up)) + length(setdiff(detected_down_edgeR, meta_down))
true_negatives <- nrow(meta_data) - (length(meta_up) + length(meta_down)) - false_positives
false_negatives <- length(setdiff(meta_up, detected_up_edgeR)) + length(setdiff(meta_down, detected_down_edgeR))
accuracy <- (true_positives + true_negatives) / nrow(meta_data)
precision <- true_positives / (true_positives + false_positives)
recall <- true_positives / (true_positives + false_negatives)
f1_score <- 2 * ((precision * recall) / (precision + recall))
FDR <- false_positives/(true_positives+false_positives) ## Precision = 1- FDR 

metrics_df_8 <- data.frame(
  Threshold = p_value,
  True_Positives = true_positives,
  False_Positives = false_positives,
  True_Negatives = true_negatives,
  False_Negatives = false_negatives,
  Accuracy = accuracy,
  Precision = precision,
  Recall = recall,
  F1_Score = f1_score,
  FDR = FDR,
  Experiment = "9_1000_0"
)

# Combine the existing data frame and the new data frame
combined_data <- rbind(metrics_df_1, metrics_df_2,metrics_df_3,metrics_df_4,metrics_df_5,metrics_df_6,metrics_df_7,metrics_df_8)
combined_data
# Write the combined data frame back to the CSV file with append = TRUE
write.csv(combined_data, "Working Directory/Output/Threshold_Edger_9_1000_0.csv", row.names = FALSE)

