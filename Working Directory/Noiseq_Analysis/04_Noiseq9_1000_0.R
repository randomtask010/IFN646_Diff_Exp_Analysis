# Step 1: Call Library (NOISeq)
library(NOISeq)

# Step 2: Load Dataset
count_data <- read.table("synthetic_data/9_1000_0.tsv", header=TRUE, row.names=1)

# Step 3: Create a conditions factor and a matrix of data
condition <- rep(c("condition1", "condition2"), each = 9)
sample_names <- paste(condition, 1:18, sep="_") # This will create unique names like "condition1_1", "condition1_2", etc.
colnames(count_data) <- sample_names

# Creating a data frame for the sample conditions
sampleinfo <- data.frame(row.names=sample_names, condition=condition)
sampleinfo

# Step 4: Create NOISeq object
mydata <- readData(data=count_data, factors=sampleinfo)
mydata

# Step 5: Differential expression analysis with NOISeq
results <- noiseq(mydata, k = 0.5, norm = "tmm", replicates = "technical", factor = "condition")


# Step 6: Extract differentially expressed genes(total differentially exprssed genes)
DEgenes <- degenes(results, q = 0.8, M = NULL)  
DEgenes

# Step 7: Load Metadata
meta_data <- read.table("synthetic_data/9_1000_0_meta.tsv", header=TRUE, row.names=1)
meta_data

# Step 8: Merge with Meta Data
annotated_results <- merge(as.data.frame(DEgenes), meta_data, by="row.names", all.x=TRUE)
annotated_results
rownames(annotated_results) <- annotated_results$Row.names
annotated_results$Row.names <- NULL


# Step 9: Comparison
#A gene is upregulated in condition 1 compared to condition 2 (equals when M = "up"), 
#it is the same thing as saying it is downregulated in condition 2 compared to condition 1.

# Downregulated in condition 2 compared to condition 1
detected_down = degenes(results, q = 0.8, M = "up")
detected_down

annotated_results_down <- merge(as.data.frame(detected_down), meta_data, by="row.names", all.x=TRUE)
annotated_results_down
rownames(annotated_results_down) <- annotated_results_down$Row.names
annotated_results_down$Row.names <- NULL
detected_down_NOISeq <- rownames(annotated_results_down)
detected_down_NOISeq

#from differentially expressed genes, actually downregulated according to metadata
meta_down <- rownames(annotated_results[annotated_results$downregulation == 1,])
meta_down

common_down <- intersect(detected_down_NOISeq, meta_down)
length(common_down)

#Upregulated
#A gene is downregulated in condition 1 compared to condition 2 (equals when M = "down"), 
#it is the same thing as saying it is upregulated in condition 2 compared to condition 1.

# upregulated in condition 2 compared to condition 1
detected_up = degenes(results, q = 0.8, M = "down")
detected_up

annotated_results_up <- merge(as.data.frame(detected_up), meta_data, by="row.names", all.x=TRUE)
annotated_results_up
rownames(annotated_results_up) <- annotated_results_up$Row.names
annotated_results_up$Row.names <- NULL
detected_up_NOISeq <- rownames(annotated_results_up)
detected_up_NOISeq

#from differentially expressed genes, actually downregulated according to metadata
meta_up <- rownames(annotated_results[annotated_results$upregulation == 1,])
meta_up

common_up <- intersect(detected_up_NOISeq, meta_up)
length(common_up)

# Step 10: Summarize outliers
outliers_up <- setdiff(detected_up_NOISeq, meta_up)
outliers_up
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

FDR <- false_positives/nrow(DEgenes)

nrow(DEgenes) #total differentilly expressed genes identified by noiseq
nrow(detected_down) #total down regulated differentilly expressed genes identified by noiseq
nrow(detected_up) #total up regulated differentilly expressed genes identified by noiseq
true_positives
false_positives
true_negatives
false_negatives
accuracy
precision
recall
f1_score
FDR

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

