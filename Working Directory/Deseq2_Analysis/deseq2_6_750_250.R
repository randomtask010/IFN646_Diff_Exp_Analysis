# Step 1: Call Library (DESeq2)
library(DESeq2)

# Step 2: Load Dataset
count_data <- read.table("RAW data/6_750_250.tsv", header=TRUE, row.names=1)
sample_info <- data.frame(groups = factor(rep(1:2, each=6)))

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

res <- res[!is.na(res$padj < 0.05), ]

# Step 8: Print top genes
head(res)

# Step 11: Load Metadata
meta_data <- read.table("RAW data/6_750_250_meta.tsv", header=TRUE, row.names=1)


# Extract DESeq2 results
DESEQ2_results <- res

# Merge DESeq2 results with metadata
annotated_results <- merge(as.data.frame(DESEQ2_results), meta_data, by="row.names", all.x=TRUE)
rownames(annotated_results) <- annotated_results$Row.names
annotated_results$Row.names <- NULL

# Step 12: Comparison using LogFC
# Upregulated
detected_up_DESeq2 <- rownames(annotated_results[annotated_results$log2FoldChange > 1 & annotated_results$padj < 0.05,])
meta_up <- rownames(annotated_results[annotated_results$upregulation == 1,])
common_up <- intersect(detected_up_DESeq2, meta_up)

# Downregulated
detected_down_DESeq2 <- rownames(annotated_results[annotated_results$log2FoldChange < -1 & annotated_results$padj < 0.05,])
meta_down <- rownames(annotated_results[annotated_results$downregulation == 1,])
common_down <- intersect(detected_down_DESeq2, meta_down)

# Step 13: Summarize outliers
outliers_up <- setdiff(detected_up_DESeq2, meta_up)
write.csv(outliers_up, "Working Directory/Output/deseq2_6_750_250_outliers_upregulated.csv", row.names = FALSE)

outliers_down <- setdiff(detected_down_DESeq2, meta_down)
write.csv(outliers_down, "Working Directory/Output/deseq2_6_750_250_outliers_downregulated.csv", row.names = FALSE)

# Step 14: Accuracy and Precision Matrix
true_positives <- length(common_up) + length(common_down)
false_positives <- length(setdiff(detected_up_DESeq2, meta_up)) + length(setdiff(detected_down_DESeq2, meta_down))
true_negatives <- nrow(meta_data) - (length(meta_up) + length(meta_down)) - false_positives
false_negatives <- length(setdiff(meta_up, detected_up_DESeq2)) + length(setdiff(meta_down, detected_down_DESeq2))
accuracy <- (true_positives + true_negatives) / nrow(meta_data)
precision <- true_positives / (true_positives + false_positives)
recall <- true_positives / (true_positives + false_negatives)
f1_score <- 2 * ((precision * recall) / (precision + recall))

# Step 15: Output Metrics to CSV for further analysis outside R
metrics_df <- data.frame(
  True_Positives = true_positives,
  False_Positives = false_positives,
  True_Negatives = true_negatives,
  False_Negatives = false_negatives,
  Accuracy = accuracy,
  Precision = precision,
  Recall = recall,
  F1_Score = f1_score
)
write.csv(metrics_df, "Working Directory/Output/Metrics_deseq2_6_750_250.csv", row.names = FALSE)