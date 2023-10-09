# Step 1: Call Library (NOISeq)
library(NOISeq)

# Step 2: Load Dataset
count_data <- read.table("9_1000_0.tsv", header=TRUE, row.names=1)

# Step 3: Create a conditions factor and a matrix of data
condition <- rep(c("condition1", "condition2"), each = 9)
sample_names <- paste(condition, 1:18, sep="_") # This will create unique names like "condition1_1", "condition1_2", etc.
colnames(count_data) <- sample_names

# Creating a data frame for the sample conditions
sampleinfo <- data.frame(row.names=sample_names, condition=condition)

# Step 4: Create NOISeq object
mydata <- readData(data=count_data, factors=sampleinfo)

# Step 5: Differential expression analysis with NOISeq
results <- noiseq(mydata, k = 0.5, norm = "n", replicates = "technical", factor = "condition")


# Step 6: Extract differentially expressed genes
DEgenes <- degenes(results, q = 0.9, M = NULL)  # 90% probability by default


# Step 7: Load Metadata
meta_data <- read.table("9_1000_0_meta.tsv", header=TRUE, row.names=1)

# Step 8: Merge with Meta Data
annotated_results <- merge(as.data.frame(DEgenes), meta_data, by="row.names", all.x=TRUE)
rownames(annotated_results) <- annotated_results$Row.names
annotated_results$Row.names <- NULL

# Step 9: Comparison
# Upregulated
detected_up_NOISeq <- rownames(annotated_results[annotated_results$log2FC > 0,])
meta_up <- rownames(annotated_results[annotated_results$upregulation == 1,])
common_up <- intersect(detected_up_NOISeq, meta_up)

# Downregulated
detected_down_NOISeq <- rownames(annotated_results[annotated_results$log2FC < 0,])
meta_down <- rownames(annotated_results[annotated_results$downregulation == 1,])
common_down <- intersect(detected_down_NOISeq, meta_down)

# Step 10: Summarize outliers
outliers_up <- setdiff(detected_up_NOISeq, meta_up)
write.csv(outliers_up, "NOISeq2_9_1000_0_outliers_upregulated.csv", row.names = FALSE)
outliers_down <- setdiff(detected_down_NOISeq, meta_down)
write.csv(outliers_down, "NOISeq2_9_1000_0_outliers_downregulated.csv", row.names = FALSE)

# Step 11: Accuracy and Precision Matrix
true_positives <- length(common_up) + length(common_down)
false_positives <- length(setdiff(detected_up_NOISeq, meta_up)) + length(setdiff(detected_down_NOISeq, meta_down))
true_negatives <- nrow(meta_data) - (length(meta_up) + length(meta_down)) - false_positives
false_negatives <- length(setdiff(meta_up, detected_up_NOISeq)) + length(setdiff(meta_down, detected_down_NOISeq))
accuracy <- (true_positives + true_negatives) / nrow(meta_data)
precision <- true_positives / (true_positives + false_positives)
recall <- true_positives / (true_positives + false_negatives)
f1_score <- 2 * ((precision * recall) / (precision + recall))

# Step 12: Venn Diagram Creation
library(VennDiagram)

# Venn diagram for upregulated genes
venn.diagram(
  x = list(NOISeq = detected_up_NOISeq, meta = meta_up),
  category.names = c("NOISeq detected up", "Metadata up"),
  output = TRUE,
  filename = "NOISeq2_9_1000_0_venn_upregulated.png",
  output.type = "png",
  imagetype = "png",
  resolution = 300
)

# Venn diagram for downregulated genes
if(length(detected_down_NOISeq) > 0 & length(meta_down) > 0) {
  venn.diagram(
    x = list(NOISeq = detected_down_NOISeq, meta = meta_down),
    category.names = c("NOISeq detected down", "Metadata down"),
    output = TRUE,
    filename = "NOISeq2_9_1000_0_venn_downregulated.png",
    output.type = "png",
    imagetype = "png",
    resolution = 300
  )
} else {
  cat("Skipping Venn diagram for downregulated genes due to empty sets.\n")
}

