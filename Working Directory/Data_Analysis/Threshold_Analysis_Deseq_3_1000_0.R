# Step 1: Call Library (DESeq2)
library(DESeq2)

# Step 2: Load Dataset
count_data <- read.table("RAW data/3_1000_0.tsv", header=TRUE, row.names=1)
sample_info <- data.frame(groups = factor(rep(1:2, each=3)))

# Step 3: Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = sample_info,
                              design = ~ groups)

# Step 4: Filter low count genes (optional)
dds <- dds[ rowSums(counts(dds)) >= 10, ]

#normalization 
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)

# Step 5: Perform differential expression analysis
dds <- DESeq(dds)

# Step 6: Extract results
res <- results(dds)

# Step 7: Summarize and filter results (optional)
summary(res)
# Subset results with padj < 0.05, while excluding NAs
padj_val <- 0.05
res <- res[!is.na(res$padj < padj_val), ]

# Step 8: Print top genes
head(res)

# Step 11: Load Metadata
meta_data <- read.table("RAW data/3_1000_0_meta.tsv", header=TRUE, row.names=1)


# Extract DESeq2 results
DESEQ2_results <- res

# Merge DESeq2 results with metadata
annotated_results <- merge(as.data.frame(DESEQ2_results), meta_data, by="row.names", all.x=TRUE)
rownames(annotated_results) <- annotated_results$Row.names
annotated_results$Row.names <- NULL

# Step 12: Comparison using LogFC
# Upregulated
detected_up_DESeq2 <- rownames(annotated_results[annotated_results$log2FoldChange > 1 & annotated_results$padj < padj_val,])
meta_up <- rownames(annotated_results[annotated_results$upregulation == 1,])
common_up <- intersect(detected_up_DESeq2, meta_up)

# Downregulated
detected_down_DESeq2 <- rownames(annotated_results[annotated_results$log2FoldChange < -1 & annotated_results$padj < padj_val,])
meta_down <- rownames(annotated_results[annotated_results$downregulation == 1,])
common_down <- intersect(detected_down_DESeq2, meta_down)

# Step 14: Accuracy and Precision Matrix
true_positives <- length(common_up) + length(common_down)
false_positives <- length(setdiff(detected_up_DESeq2, meta_up)) + length(setdiff(detected_down_DESeq2, meta_down))
true_negatives <- nrow(meta_data) - (length(meta_up) + length(meta_down)) - false_positives
false_negatives <- length(setdiff(meta_up, detected_up_DESeq2)) + length(setdiff(meta_down, detected_down_DESeq2))
accuracy <- (true_positives + true_negatives) / nrow(meta_data)
precision <- true_positives / (true_positives + false_positives)
recall <- true_positives / (true_positives + false_negatives)
f1_score <- 2 * ((precision * recall) / (precision + recall))
FDR <- false_positives/(true_positives+false_positives) ## Precision = 1- FDR 

# Step 15: Output Metrics to CSV for further analysis outside R
metrics_df_1 <- data.frame(
  Threshold = padj_val,
  True_Positives = true_positives,
  False_Positives = false_positives,
  True_Negatives = true_negatives,
  False_Negatives = false_negatives,
  Accuracy = accuracy,
  Precision = precision,
  Recall = recall,
  F1_Score = f1_score,
  FDR = FDR,
  Experiment = "3_1000_0"
)

# Subset results with padj < 0.04, while excluding NAs
padj_val <- 0.04
res <- res[!is.na(res$padj < padj_val), ]

DESEQ2_results <- res

annotated_results <- merge(as.data.frame(DESEQ2_results), meta_data, by="row.names", all.x=TRUE)
rownames(annotated_results) <- annotated_results$Row.names
annotated_results$Row.names <- NULL

# Upregulated
detected_up_DESeq2 <- rownames(annotated_results[annotated_results$log2FoldChange > 1 & annotated_results$padj < padj_val,])
meta_up <- rownames(annotated_results[annotated_results$upregulation == 1,])
common_up <- intersect(detected_up_DESeq2, meta_up)

# Downregulated
detected_down_DESeq2 <- rownames(annotated_results[annotated_results$log2FoldChange < -1 & annotated_results$padj < padj_val,])
meta_down <- rownames(annotated_results[annotated_results$downregulation == 1,])
common_down <- intersect(detected_down_DESeq2, meta_down)

true_positives <- length(common_up) + length(common_down)
false_positives <- length(setdiff(detected_up_DESeq2, meta_up)) + length(setdiff(detected_down_DESeq2, meta_down))
true_negatives <- nrow(meta_data) - (length(meta_up) + length(meta_down)) - false_positives
false_negatives <- length(setdiff(meta_up, detected_up_DESeq2)) + length(setdiff(meta_down, detected_down_DESeq2))
accuracy <- (true_positives + true_negatives) / nrow(meta_data)
precision <- true_positives / (true_positives + false_positives)
recall <- true_positives / (true_positives + false_negatives)
f1_score <- 2 * ((precision * recall) / (precision + recall))
FDR <- false_positives/(true_positives+false_positives) ## Precision = 1- FDR 

metrics_df_2 <- data.frame(
  Threshold = padj_val,
  True_Positives = true_positives,
  False_Positives = false_positives,
  True_Negatives = true_negatives,
  False_Negatives = false_negatives,
  Accuracy = accuracy,
  Precision = precision,
  Recall = recall,
  F1_Score = f1_score,
  FDR = FDR,
  Experiment = "3_1000_0"
)

# Subset results with padj < 0.03, while excluding NAs
padj_val <- 0.03
res <- res[!is.na(res$padj < padj_val), ]

DESEQ2_results <- res

annotated_results <- merge(as.data.frame(DESEQ2_results), meta_data, by="row.names", all.x=TRUE)
rownames(annotated_results) <- annotated_results$Row.names
annotated_results$Row.names <- NULL

# Upregulated
detected_up_DESeq2 <- rownames(annotated_results[annotated_results$log2FoldChange > 1 & annotated_results$padj < padj_val,])
meta_up <- rownames(annotated_results[annotated_results$upregulation == 1,])
common_up <- intersect(detected_up_DESeq2, meta_up)

# Downregulated
detected_down_DESeq2 <- rownames(annotated_results[annotated_results$log2FoldChange < -1 & annotated_results$padj < padj_val,])
meta_down <- rownames(annotated_results[annotated_results$downregulation == 1,])
common_down <- intersect(detected_down_DESeq2, meta_down)

true_positives <- length(common_up) + length(common_down)
false_positives <- length(setdiff(detected_up_DESeq2, meta_up)) + length(setdiff(detected_down_DESeq2, meta_down))
true_negatives <- nrow(meta_data) - (length(meta_up) + length(meta_down)) - false_positives
false_negatives <- length(setdiff(meta_up, detected_up_DESeq2)) + length(setdiff(meta_down, detected_down_DESeq2))
accuracy <- (true_positives + true_negatives) / nrow(meta_data)
precision <- true_positives / (true_positives + false_positives)
recall <- true_positives / (true_positives + false_negatives)
f1_score <- 2 * ((precision * recall) / (precision + recall))
FDR <- false_positives/(true_positives+false_positives) ## Precision = 1- FDR 

metrics_df_3 <- data.frame(
  Threshold = padj_val,
  True_Positives = true_positives,
  False_Positives = false_positives,
  True_Negatives = true_negatives,
  False_Negatives = false_negatives,
  Accuracy = accuracy,
  Precision = precision,
  Recall = recall,
  F1_Score = f1_score,
  FDR = FDR,
  Experiment = "3_1000_0"
)

# Subset results with padj < 0.06, while excluding NAs
padj_val <- 0.06
res <- res[!is.na(res$padj < padj_val), ]

DESEQ2_results <- res

annotated_results <- merge(as.data.frame(DESEQ2_results), meta_data, by="row.names", all.x=TRUE)
rownames(annotated_results) <- annotated_results$Row.names
annotated_results$Row.names <- NULL

# Upregulated
detected_up_DESeq2 <- rownames(annotated_results[annotated_results$log2FoldChange > 1 & annotated_results$padj < padj_val,])
meta_up <- rownames(annotated_results[annotated_results$upregulation == 1,])
common_up <- intersect(detected_up_DESeq2, meta_up)

# Downregulated
detected_down_DESeq2 <- rownames(annotated_results[annotated_results$log2FoldChange < -1 & annotated_results$padj < padj_val,])
meta_down <- rownames(annotated_results[annotated_results$downregulation == 1,])
common_down <- intersect(detected_down_DESeq2, meta_down)

true_positives <- length(common_up) + length(common_down)
false_positives <- length(setdiff(detected_up_DESeq2, meta_up)) + length(setdiff(detected_down_DESeq2, meta_down))
true_negatives <- nrow(meta_data) - (length(meta_up) + length(meta_down)) - false_positives
false_negatives <- length(setdiff(meta_up, detected_up_DESeq2)) + length(setdiff(meta_down, detected_down_DESeq2))
accuracy <- (true_positives + true_negatives) / nrow(meta_data)
precision <- true_positives / (true_positives + false_positives)
recall <- true_positives / (true_positives + false_negatives)
f1_score <- 2 * ((precision * recall) / (precision + recall))
FDR <- false_positives/(true_positives+false_positives) ## Precision = 1- FDR 

metrics_df_4 <- data.frame(
  Threshold = padj_val,
  True_Positives = true_positives,
  False_Positives = false_positives,
  True_Negatives = true_negatives,
  False_Negatives = false_negatives,
  Accuracy = accuracy,
  Precision = precision,
  Recall = recall,
  F1_Score = f1_score,
  FDR = FDR,
  Experiment = "3_1000_0"
)

# Subset results with padj < 0.07, while excluding NAs
padj_val <- 0.07
res <- res[!is.na(res$padj < padj_val), ]

DESEQ2_results <- res

annotated_results <- merge(as.data.frame(DESEQ2_results), meta_data, by="row.names", all.x=TRUE)
rownames(annotated_results) <- annotated_results$Row.names
annotated_results$Row.names <- NULL

# Upregulated
detected_up_DESeq2 <- rownames(annotated_results[annotated_results$log2FoldChange > 1 & annotated_results$padj < padj_val,])
meta_up <- rownames(annotated_results[annotated_results$upregulation == 1,])
common_up <- intersect(detected_up_DESeq2, meta_up)

# Downregulated
detected_down_DESeq2 <- rownames(annotated_results[annotated_results$log2FoldChange < -1 & annotated_results$padj < padj_val,])
meta_down <- rownames(annotated_results[annotated_results$downregulation == 1,])
common_down <- intersect(detected_down_DESeq2, meta_down)

true_positives <- length(common_up) + length(common_down)
false_positives <- length(setdiff(detected_up_DESeq2, meta_up)) + length(setdiff(detected_down_DESeq2, meta_down))
true_negatives <- nrow(meta_data) - (length(meta_up) + length(meta_down)) - false_positives
false_negatives <- length(setdiff(meta_up, detected_up_DESeq2)) + length(setdiff(meta_down, detected_down_DESeq2))
accuracy <- (true_positives + true_negatives) / nrow(meta_data)
precision <- true_positives / (true_positives + false_positives)
recall <- true_positives / (true_positives + false_negatives)
f1_score <- 2 * ((precision * recall) / (precision + recall))
FDR <- false_positives/(true_positives+false_positives) ## Precision = 1- FDR 

metrics_df_5 <- data.frame(
  Threshold = padj_val,
  True_Positives = true_positives,
  False_Positives = false_positives,
  True_Negatives = true_negatives,
  False_Negatives = false_negatives,
  Accuracy = accuracy,
  Precision = precision,
  Recall = recall,
  F1_Score = f1_score,
  FDR = FDR,
  Experiment = "3_1000_0"
)

# Subset results with padj < 0.08, while excluding NAs
padj_val <- 0.08
res <- res[!is.na(res$padj < padj_val), ]

DESEQ2_results <- res

annotated_results <- merge(as.data.frame(DESEQ2_results), meta_data, by="row.names", all.x=TRUE)
rownames(annotated_results) <- annotated_results$Row.names
annotated_results$Row.names <- NULL

# Upregulated
detected_up_DESeq2 <- rownames(annotated_results[annotated_results$log2FoldChange > 1 & annotated_results$padj < padj_val,])
meta_up <- rownames(annotated_results[annotated_results$upregulation == 1,])
common_up <- intersect(detected_up_DESeq2, meta_up)

# Downregulated
detected_down_DESeq2 <- rownames(annotated_results[annotated_results$log2FoldChange < -1 & annotated_results$padj < padj_val,])
meta_down <- rownames(annotated_results[annotated_results$downregulation == 1,])
common_down <- intersect(detected_down_DESeq2, meta_down)

true_positives <- length(common_up) + length(common_down)
false_positives <- length(setdiff(detected_up_DESeq2, meta_up)) + length(setdiff(detected_down_DESeq2, meta_down))
true_negatives <- nrow(meta_data) - (length(meta_up) + length(meta_down)) - false_positives
false_negatives <- length(setdiff(meta_up, detected_up_DESeq2)) + length(setdiff(meta_down, detected_down_DESeq2))
accuracy <- (true_positives + true_negatives) / nrow(meta_data)
precision <- true_positives / (true_positives + false_positives)
recall <- true_positives / (true_positives + false_negatives)
f1_score <- 2 * ((precision * recall) / (precision + recall))
FDR <- false_positives/(true_positives+false_positives) ## Precision = 1- FDR 

metrics_df_6 <- data.frame(
  Threshold = padj_val,
  True_Positives = true_positives,
  False_Positives = false_positives,
  True_Negatives = true_negatives,
  False_Negatives = false_negatives,
  Accuracy = accuracy,
  Precision = precision,
  Recall = recall,
  F1_Score = f1_score,
  FDR = FDR,
  Experiment = "3_1000_0"
)

# Subset results with padj < 0.09, while excluding NAs
padj_val <- 0.09
res <- res[!is.na(res$padj < padj_val), ]

DESEQ2_results <- res

annotated_results <- merge(as.data.frame(DESEQ2_results), meta_data, by="row.names", all.x=TRUE)
rownames(annotated_results) <- annotated_results$Row.names
annotated_results$Row.names <- NULL

# Upregulated
detected_up_DESeq2 <- rownames(annotated_results[annotated_results$log2FoldChange > 1 & annotated_results$padj < padj_val,])
meta_up <- rownames(annotated_results[annotated_results$upregulation == 1,])
common_up <- intersect(detected_up_DESeq2, meta_up)

# Downregulated
detected_down_DESeq2 <- rownames(annotated_results[annotated_results$log2FoldChange < -1 & annotated_results$padj < padj_val,])
meta_down <- rownames(annotated_results[annotated_results$downregulation == 1,])
common_down <- intersect(detected_down_DESeq2, meta_down)

true_positives <- length(common_up) + length(common_down)
false_positives <- length(setdiff(detected_up_DESeq2, meta_up)) + length(setdiff(detected_down_DESeq2, meta_down))
true_negatives <- nrow(meta_data) - (length(meta_up) + length(meta_down)) - false_positives
false_negatives <- length(setdiff(meta_up, detected_up_DESeq2)) + length(setdiff(meta_down, detected_down_DESeq2))
accuracy <- (true_positives + true_negatives) / nrow(meta_data)
precision <- true_positives / (true_positives + false_positives)
recall <- true_positives / (true_positives + false_negatives)
f1_score <- 2 * ((precision * recall) / (precision + recall))
FDR <- false_positives/(true_positives+false_positives) ## Precision = 1- FDR 

metrics_df_7 <- data.frame(
  Threshold = padj_val,
  True_Positives = true_positives,
  False_Positives = false_positives,
  True_Negatives = true_negatives,
  False_Negatives = false_negatives,
  Accuracy = accuracy,
  Precision = precision,
  Recall = recall,
  F1_Score = f1_score,
  FDR = FDR,
  Experiment = "3_1000_0"
)

# Subset results with padj < 0.09, while excluding NAs
padj_val <- 0.02
res <- res[!is.na(res$padj < padj_val), ]

DESEQ2_results <- res

annotated_results <- merge(as.data.frame(DESEQ2_results), meta_data, by="row.names", all.x=TRUE)
rownames(annotated_results) <- annotated_results$Row.names
annotated_results$Row.names <- NULL

# Upregulated
detected_up_DESeq2 <- rownames(annotated_results[annotated_results$log2FoldChange > 1 & annotated_results$padj < padj_val,])
meta_up <- rownames(annotated_results[annotated_results$upregulation == 1,])
common_up <- intersect(detected_up_DESeq2, meta_up)

# Downregulated
detected_down_DESeq2 <- rownames(annotated_results[annotated_results$log2FoldChange < -1 & annotated_results$padj < padj_val,])
meta_down <- rownames(annotated_results[annotated_results$downregulation == 1,])
common_down <- intersect(detected_down_DESeq2, meta_down)

true_positives <- length(common_up) + length(common_down)
false_positives <- length(setdiff(detected_up_DESeq2, meta_up)) + length(setdiff(detected_down_DESeq2, meta_down))
true_negatives <- nrow(meta_data) - (length(meta_up) + length(meta_down)) - false_positives
false_negatives <- length(setdiff(meta_up, detected_up_DESeq2)) + length(setdiff(meta_down, detected_down_DESeq2))
accuracy <- (true_positives + true_negatives) / nrow(meta_data)
precision <- true_positives / (true_positives + false_positives)
recall <- true_positives / (true_positives + false_negatives)
f1_score <- 2 * ((precision * recall) / (precision + recall))
FDR <- false_positives/(true_positives+false_positives) ## Precision = 1- FDR 

metrics_df_8 <- data.frame(
  Threshold = padj_val,
  True_Positives = true_positives,
  False_Positives = false_positives,
  True_Negatives = true_negatives,
  False_Negatives = false_negatives,
  Accuracy = accuracy,
  Precision = precision,
  Recall = recall,
  F1_Score = f1_score,
  FDR = FDR,
  Experiment = "3_1000_0"
)

# Combine the existing data frame and the new data frame
combined_data <- rbind(metrics_df_1, metrics_df_2,metrics_df_3,metrics_df_4,metrics_df_5,metrics_df_6,metrics_df_7,metrics_df_8)
combined_data
# Write the combined data frame back to the CSV file with append = TRUE
write.csv(combined_data, "Working Directory/Output/Threshold_deseq_3_1000_0.csv", row.names = FALSE)

