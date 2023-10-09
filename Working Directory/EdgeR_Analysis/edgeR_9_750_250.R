# Step 1: Call Library (edgeR)
library(edgeR)

# Step 2: Load Dataset
count_data <- read.table("RAW data/9_750_250.tsv", header=TRUE, row.names=1)
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
decided <- decideTestsDGE(et, adjust.method="BH", p.value=0.05)
summary(decided)

# Step 10: Print top genes
topTags(et)

# Step 11: Load Metadata
meta_data <- read.table("RAW data/9_750_250_meta.tsv", header=TRUE, row.names=1)
edgeR_results <- topTags(et, n=Inf)
annotated_results <- merge(edgeR_results, meta_data, by="row.names", all.x=TRUE)
rownames(annotated_results) <- annotated_results$Row.names
annotated_results$Row.names <- NULL

# Step 12: Comparison using LogFC
# Upregulated
detected_up_edgeR <- rownames(annotated_results[annotated_results$logFC > 0 & annotated_results$FDR < 0.05,])
meta_up <- rownames(annotated_results[annotated_results$upregulation == 1,])
common_up <- intersect(detected_up_edgeR, meta_up)
# Downregulated
detected_down_edgeR <- rownames(annotated_results[annotated_results$logFC < 0 & annotated_results$FDR < 0.05,])
meta_down <- rownames(annotated_results[annotated_results$downregulation == 1,])
common_down <- intersect(detected_down_edgeR, meta_down)

# Step 13: Summarize outliers
outliers_up <- setdiff(detected_up_edgeR, meta_up)
write.csv(outliers_up, "Working Directory/Output/edgeR_9_750_250_outliers_upregulated.csv", row.names = FALSE)
outliers_down <- setdiff(detected_down_edgeR, meta_down)
write.csv(outliers_down, "Working Directory/Output/edgeR_9_750_250_outliers_downregulated.csv", row.names = FALSE)

# Step 14: Accuracy and Precision Matrix
true_positives <- length(common_up) + length(common_down)
false_positives <- length(setdiff(detected_up_edgeR, meta_up)) + length(setdiff(detected_down_edgeR, meta_down))
true_negatives <- nrow(meta_data) - (length(meta_up) + length(meta_down)) - false_positives
false_negatives <- length(setdiff(meta_up, detected_up_edgeR)) + length(setdiff(meta_down, detected_down_edgeR))
accuracy <- (true_positives + true_negatives) / nrow(meta_data)
precision <- true_positives / (true_positives + false_positives)
recall <- true_positives / (true_positives + false_negatives)
f1_score <- 2 * ((precision * recall) / (precision + recall))
