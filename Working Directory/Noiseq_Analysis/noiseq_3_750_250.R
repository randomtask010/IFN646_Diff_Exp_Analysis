# Step 1: Call Library (NOISeq)
library(NOISeq)

# Step 2: Load Dataset
count_data <- read.table("RAW data/3_750_250.tsv", header=TRUE, row.names=1)

# Step 3: Create a conditions factor and a matrix of data
condition <- rep(c("condition1", "condition2"), each = 3)
sample_names <- paste(condition, 1:6, sep="_") # This will create unique names like "condition1_1", "condition1_2", etc.
colnames(count_data) <- sample_names

# Creating a data frame for the sample conditions
sampleinfo <- data.frame(row.names=sample_names, condition=condition)
sampleinfo

#Filtering data
#Excluding features with low counts improves differential expression results, since noise in the data is reduced
#counts per million (CPM) method is used to filtering
count_data = filtered.data(count_data, factor = sampleinfo$condition, norm=FALSE, depth=NULL, method=1, cv.cutoff=100, cpm=1, p.adj="fdr")
count_data

# Step 4: Create NOISeq object
mydata <- readData(data=count_data, factors=sampleinfo)
mydata

# Step 5: Differential expression analysis with NOISeq
results <- noiseq(mydata, k = 0.5, norm = "tmm", replicates = "technical", factor = "condition")


# Step 6: Extract differentially expressed genes(total differentially exprssed genes)
DEgenes <- degenes(results, q = 0.8, M = NULL)  
DEgenes

# Step 7: Load Metadata
meta_data <- read.table("RAW data/3_750_250_meta.tsv", header=TRUE, row.names=1)
meta_data

# Step 8: Merge with Meta Data
annotated_results <- merge(as.data.frame(DEgenes), meta_data, by="row.names", all.x=TRUE)
annotated_results
rownames(annotated_results) <- annotated_results$Row.names
annotated_results$Row.names <- NULL


# Step 9: Comparison
#A gene is upregulated in condition 1 compared to condition 2 (equals when M = "up"), 
#it is the same thing as saying it is downregulated in condition 2 compared to condition 1.
#detected_down represent the down regulated genes in condition 2 compared to condition 1.
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
#detected_up represent the up regulated genes in condition 2 compared to condition 1.
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
write.csv(outliers_up, "Working Directory/Output/noiseq_3_750_250_outliers_upregulated.csv", row.names = FALSE)
outliers_down <- setdiff(detected_down_NOISeq, meta_down)
write.csv(outliers_down, "Working Directory/Output/noiseq_3_750_250_outliers_downregulated.csv", row.names = FALSE)

# Step 11: Accuracy and Precision Matrix
true_positives <- length(common_up) + length(common_down)
false_positives <- length(setdiff(detected_up_NOISeq, meta_up)) + length(setdiff(detected_down_NOISeq, meta_down))
true_negatives <- nrow(meta_data) - (length(meta_up) + length(meta_down)) - false_positives
false_negatives <- length(setdiff(meta_up, detected_up_NOISeq)) + length(setdiff(meta_down, detected_down_NOISeq))
accuracy <- (true_positives + true_negatives) / nrow(meta_data)
precision <- true_positives / (true_positives + false_positives)
recall <- true_positives / (true_positives + false_negatives)
f1_score <- 2 * ((precision * recall) / (precision + recall))

# Step 12: Output Metrics to CSV for further analysis outside R
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
write.csv(metrics_df, "Working Directory/Output/Metrics_noiseq_3_750_250.csv", row.names = FALSE)

# BELOW TO BE Included in overall metric output if needed or deleted if used for debugging - NISH to action
# could be useful when expanding metrics to include time for code to execute - analysis output as a proxy for efficency 

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

