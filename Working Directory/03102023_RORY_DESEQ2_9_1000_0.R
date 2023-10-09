# Step 1: Call Library (DESeq2)
library(DESeq2)

# Step 2: Load Dataset
count_data <- read.table("9_1000_0.tsv", header=TRUE, row.names=1)
col_data <- data.frame(groups = factor(rep(1:2, each=9)))

# Step 3: Create DESeqDataSet Object
dds <- DESeqDataSetFromMatrix(countData = count_data, colData = col_data, design = ~groups)

# Step 4: Pre-filtering to remove low count genes
keep <- rowSums(counts(dds)) >= 10  # for instance, you can adjust the threshold as needed
dds <- dds[keep,]

# Step 5: Differential Expression Analysis
dds <- DESeq(dds)

# Step 6: Extract Differential Expression data
res <- results(dds)

# Step 7: Order by adjusted p-value and extract top genes
res_ordered <- res[order(res$padj), ]

# Step 8: Load Metadata
meta_data <- read.table("9_1000_0_meta.tsv", header=TRUE, row.names=1)
annotated_results <- merge(as.data.frame(res_ordered), meta_data, by="row.names", all.x=TRUE)
rownames(annotated_results) <- annotated_results$Row.names
annotated_results$Row.names <- NULL

# Step 9: Comparison using Log2FoldChange
# Upregulated
detected_up_DESeq2 <- rownames(annotated_results[annotated_results$log2FoldChange > 0 & annotated_results$padj < 0.05,])
meta_up <- rownames(annotated_results[annotated_results$upregulation == 1,])
common_up <- intersect(detected_up_DESeq2, meta_up)

# Downregulated
detected_down_DESeq2 <- rownames(annotated_results[annotated_results$log2FoldChange < 0 & annotated_results$padj < 0.05,])
meta_down <- rownames(annotated_results[annotated_results$downregulation == 1,])
common_down <- intersect(detected_down_DESeq2, meta_down)

# Step 10: Summarize outliers
outliers_up <- setdiff(detected_up_DESeq2, meta_up)
write.csv(outliers_up, "DESEQ2_9_1000_0_outliers_upregulated.csv", row.names = FALSE)
outliers_down <- setdiff(detected_down_DESeq2, meta_down)
write.csv(outliers_down, "DESEQ2_9_1000_0_outliers_downregulated.csv", row.names = FALSE)

# Step 11: Accuracy and Precision Matrix
true_positives <- length(common_up) + length(common_down)
false_positives <- length(setdiff(detected_up_DESeq2, meta_up)) + length(setdiff(detected_down_DESeq2, meta_down))
true_negatives <- nrow(meta_data) - (length(meta_up) + length(meta_down)) - false_positives
false_negatives <- length(setdiff(meta_up, detected_up_DESeq2)) + length(setdiff(meta_down, detected_down_DESeq2))
accuracy <- (true_positives + true_negatives) / nrow(meta_data)
precision <- true_positives / (true_positives + false_positives)
recall <- true_positives / (true_positives + false_negatives)
f1_score <- 2 * ((precision * recall) / (precision + recall))

# Step 12: Venn Diagram Creation
library(VennDiagram)

# Venn diagram for upregulated genes
venn.diagram(
  x = list(DESeq2 = detected_up_DESeq2, meta = meta_up),
  category.names = c("DESeq2 detected up", "Metadata up"),
  output = TRUE,
  filename = "DESEQ2_9_1000_0_venn_upregulated.png",
  output.type = "png",
  imagetype = "png",
  resolution = 300
)

# Venn diagram for downregulated genes
venn.diagram(
  x = list(DESeq2 = detected_down_DESeq2, meta = meta_down),
  category.names = c("DESeq2 detected down", "Metadata down"),
  output = TRUE,
  filename = "DESEQ2_9_1000_0_venn_downregulated.png",
  output.type = "png",
  imagetype = "png",
  resolution = 300
)
