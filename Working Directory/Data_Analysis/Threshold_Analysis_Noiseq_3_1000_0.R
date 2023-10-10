# Step 1: Call Library (NOISeq)
library(NOISeq)

# Step 2: Load Dataset
count_data <- read.table("RAW data/3_1000_0.tsv", header=TRUE, row.names=1)

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

q_value <- 0.8
DEgenes <- degenes(results, q = q_value, M = NULL)  
meta_data <- read.table("RAW data/3_1000_0_meta.tsv", header=TRUE, row.names=1)

annotated_results <- merge(as.data.frame(DEgenes), meta_data, by="row.names", all.x=TRUE)
annotated_results
rownames(annotated_results) <- annotated_results$Row.names
annotated_results$Row.names <- NULL

# Downregulated in condition 2 compared to condition 1
detected_down = degenes(results, q = q_value, M = "up")
detected_down

annotated_results_down <- merge(as.data.frame(detected_down), meta_data, by="row.names", all.x=TRUE)
annotated_results_down
rownames(annotated_results_down) <- annotated_results_down$Row.names
annotated_results_down$Row.names <- NULL
detected_down_NOISeq <- rownames(annotated_results_down)
detected_down_NOISeq

#from differentially expressed genes, actually downregulated according to metadata
meta_down <- rownames(annotated_results[annotated_results$downregulation == 1,])

common_down <- intersect(detected_down_NOISeq, meta_down)

#detected_up represent the up regulated genes in condition 2 compared to condition 1.
detected_up = degenes(results, q = q_value, M = "down")

annotated_results_up <- merge(as.data.frame(detected_up), meta_data, by="row.names", all.x=TRUE)
rownames(annotated_results_up) <- annotated_results_up$Row.names
annotated_results_up$Row.names <- NULL
detected_up_NOISeq <- rownames(annotated_results_up)


#from differentially expressed genes, actually downregulated according to metadata
meta_up <- rownames(annotated_results[annotated_results$upregulation == 1,])

common_up <- intersect(detected_up_NOISeq, meta_up)

# Step 10: Summarize outliers
outliers_up <- setdiff(detected_up_NOISeq, meta_up)

# Step 11: Accuracy and Precision Matrix
true_positives <- length(common_up) + length(common_down)
false_positives <- length(setdiff(detected_up_NOISeq, meta_up)) + length(setdiff(detected_down_NOISeq, meta_down))
true_negatives <- nrow(meta_data) - (length(meta_up) + length(meta_down)) - false_positives
false_negatives <- length(setdiff(meta_up, detected_up_NOISeq)) + length(setdiff(meta_down, detected_down_NOISeq))
accuracy <- (true_positives + true_negatives) / nrow(meta_data)
precision <- true_positives / (true_positives + false_positives)
recall <- true_positives / (true_positives + false_negatives)
f1_score <- 2 * ((precision * recall) / (precision + recall))
FDR <- false_positives/(true_positives+false_positives) ## Precision = 1- FDR 

# Step 12: Output Metrics to CSV for further analysis outside R
metrics_df_1 <- data.frame(
  Threshold = q_value,
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

## Executing the same for q=0.9
q_value <- 0.9
DEgenes <- degenes(results, q = q_value, M = NULL)  

annotated_results <- merge(as.data.frame(DEgenes), meta_data, by="row.names", all.x=TRUE)
annotated_results
rownames(annotated_results) <- annotated_results$Row.names
annotated_results$Row.names <- NULL

# Downregulated in condition 2 compared to condition 1
detected_down = degenes(results, q = q_value, M = "up")
detected_down

annotated_results_down <- merge(as.data.frame(detected_down), meta_data, by="row.names", all.x=TRUE)
annotated_results_down
rownames(annotated_results_down) <- annotated_results_down$Row.names
annotated_results_down$Row.names <- NULL
detected_down_NOISeq <- rownames(annotated_results_down)
detected_down_NOISeq

#from differentially expressed genes, actually downregulated according to metadata
meta_down <- rownames(annotated_results[annotated_results$downregulation == 1,])

common_down <- intersect(detected_down_NOISeq, meta_down)

#detected_up represent the up regulated genes in condition 2 compared to condition 1.
detected_up = degenes(results, q = q_value, M = "down")

annotated_results_up <- merge(as.data.frame(detected_up), meta_data, by="row.names", all.x=TRUE)
rownames(annotated_results_up) <- annotated_results_up$Row.names
annotated_results_up$Row.names <- NULL
detected_up_NOISeq <- rownames(annotated_results_up)


#from differentially expressed genes, actually downregulated according to metadata
meta_up <- rownames(annotated_results[annotated_results$upregulation == 1,])

common_up <- intersect(detected_up_NOISeq, meta_up)

# Step 10: Summarize outliers
outliers_up <- setdiff(detected_up_NOISeq, meta_up)

# Step 11: Accuracy and Precision Matrix
true_positives <- length(common_up) + length(common_down)
false_positives <- length(setdiff(detected_up_NOISeq, meta_up)) + length(setdiff(detected_down_NOISeq, meta_down))
true_negatives <- nrow(meta_data) - (length(meta_up) + length(meta_down)) - false_positives
false_negatives <- length(setdiff(meta_up, detected_up_NOISeq)) + length(setdiff(meta_down, detected_down_NOISeq))
accuracy <- (true_positives + true_negatives) / nrow(meta_data)
precision <- true_positives / (true_positives + false_positives)
recall <- true_positives / (true_positives + false_negatives)
f1_score <- 2 * ((precision * recall) / (precision + recall))
FDR <- false_positives/(true_positives+false_positives) ## Precision = 1- FDR 

# Step 12: Output Metrics to CSV for further analysis outside R
metrics_df_2 <- data.frame(
  Threshold = q_value,
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
## Executing the same for q=0.7
q_value <- 0.7
DEgenes <- degenes(results, q = q_value, M = NULL)  

annotated_results <- merge(as.data.frame(DEgenes), meta_data, by="row.names", all.x=TRUE)
annotated_results
rownames(annotated_results) <- annotated_results$Row.names
annotated_results$Row.names <- NULL

# Downregulated in condition 2 compared to condition 1
detected_down = degenes(results, q = q_value, M = "up")
detected_down

annotated_results_down <- merge(as.data.frame(detected_down), meta_data, by="row.names", all.x=TRUE)
annotated_results_down
rownames(annotated_results_down) <- annotated_results_down$Row.names
annotated_results_down$Row.names <- NULL
detected_down_NOISeq <- rownames(annotated_results_down)
detected_down_NOISeq

#from differentially expressed genes, actually downregulated according to metadata
meta_down <- rownames(annotated_results[annotated_results$downregulation == 1,])

common_down <- intersect(detected_down_NOISeq, meta_down)

#detected_up represent the up regulated genes in condition 2 compared to condition 1.
detected_up = degenes(results, q = q_value, M = "down")

annotated_results_up <- merge(as.data.frame(detected_up), meta_data, by="row.names", all.x=TRUE)
rownames(annotated_results_up) <- annotated_results_up$Row.names
annotated_results_up$Row.names <- NULL
detected_up_NOISeq <- rownames(annotated_results_up)


#from differentially expressed genes, actually downregulated according to metadata
meta_up <- rownames(annotated_results[annotated_results$upregulation == 1,])

common_up <- intersect(detected_up_NOISeq, meta_up)

# Step 10: Summarize outliers
outliers_up <- setdiff(detected_up_NOISeq, meta_up)

# Step 11: Accuracy and Precision Matrix
true_positives <- length(common_up) + length(common_down)
false_positives <- length(setdiff(detected_up_NOISeq, meta_up)) + length(setdiff(detected_down_NOISeq, meta_down))
true_negatives <- nrow(meta_data) - (length(meta_up) + length(meta_down)) - false_positives
false_negatives <- length(setdiff(meta_up, detected_up_NOISeq)) + length(setdiff(meta_down, detected_down_NOISeq))
accuracy <- (true_positives + true_negatives) / nrow(meta_data)
precision <- true_positives / (true_positives + false_positives)
recall <- true_positives / (true_positives + false_negatives)
f1_score <- 2 * ((precision * recall) / (precision + recall))
FDR <- false_positives/(true_positives+false_positives) ## Precision = 1- FDR 

# Step 12: Output Metrics to CSV for further analysis outside R
metrics_df_3 <- data.frame(
  Threshold = q_value,
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
## Executing the same for q=0.6
q_value <- 0.6
DEgenes <- degenes(results, q = q_value, M = NULL)  

annotated_results <- merge(as.data.frame(DEgenes), meta_data, by="row.names", all.x=TRUE)
annotated_results
rownames(annotated_results) <- annotated_results$Row.names
annotated_results$Row.names <- NULL

# Downregulated in condition 2 compared to condition 1
detected_down = degenes(results, q = q_value, M = "up")
detected_down

annotated_results_down <- merge(as.data.frame(detected_down), meta_data, by="row.names", all.x=TRUE)
annotated_results_down
rownames(annotated_results_down) <- annotated_results_down$Row.names
annotated_results_down$Row.names <- NULL
detected_down_NOISeq <- rownames(annotated_results_down)
detected_down_NOISeq

#from differentially expressed genes, actually downregulated according to metadata
meta_down <- rownames(annotated_results[annotated_results$downregulation == 1,])

common_down <- intersect(detected_down_NOISeq, meta_down)

#detected_up represent the up regulated genes in condition 2 compared to condition 1.
detected_up = degenes(results, q = q_value, M = "down")

annotated_results_up <- merge(as.data.frame(detected_up), meta_data, by="row.names", all.x=TRUE)
rownames(annotated_results_up) <- annotated_results_up$Row.names
annotated_results_up$Row.names <- NULL
detected_up_NOISeq <- rownames(annotated_results_up)


#from differentially expressed genes, actually downregulated according to metadata
meta_up <- rownames(annotated_results[annotated_results$upregulation == 1,])

common_up <- intersect(detected_up_NOISeq, meta_up)

# Step 10: Summarize outliers
outliers_up <- setdiff(detected_up_NOISeq, meta_up)

# Step 11: Accuracy and Precision Matrix
true_positives <- length(common_up) + length(common_down)
false_positives <- length(setdiff(detected_up_NOISeq, meta_up)) + length(setdiff(detected_down_NOISeq, meta_down))
true_negatives <- nrow(meta_data) - (length(meta_up) + length(meta_down)) - false_positives
false_negatives <- length(setdiff(meta_up, detected_up_NOISeq)) + length(setdiff(meta_down, detected_down_NOISeq))
accuracy <- (true_positives + true_negatives) / nrow(meta_data)
precision <- true_positives / (true_positives + false_positives)
recall <- true_positives / (true_positives + false_negatives)
f1_score <- 2 * ((precision * recall) / (precision + recall))
FDR <- false_positives/(true_positives+false_positives) ## Precision = 1- FDR 

# Step 12: Output Metrics to CSV for further analysis outside R
metrics_df_4 <- data.frame(
  Threshold = q_value,
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

## Executing the same for q=0.5
q_value <- 0.5
DEgenes <- degenes(results, q = q_value, M = NULL)  

annotated_results <- merge(as.data.frame(DEgenes), meta_data, by="row.names", all.x=TRUE)
annotated_results
rownames(annotated_results) <- annotated_results$Row.names
annotated_results$Row.names <- NULL

# Downregulated in condition 2 compared to condition 1
detected_down = degenes(results, q = q_value, M = "up")
detected_down

annotated_results_down <- merge(as.data.frame(detected_down), meta_data, by="row.names", all.x=TRUE)
annotated_results_down
rownames(annotated_results_down) <- annotated_results_down$Row.names
annotated_results_down$Row.names <- NULL
detected_down_NOISeq <- rownames(annotated_results_down)
detected_down_NOISeq

#from differentially expressed genes, actually downregulated according to metadata
meta_down <- rownames(annotated_results[annotated_results$downregulation == 1,])

common_down <- intersect(detected_down_NOISeq, meta_down)

#detected_up represent the up regulated genes in condition 2 compared to condition 1.
detected_up = degenes(results, q = q_value, M = "down")

annotated_results_up <- merge(as.data.frame(detected_up), meta_data, by="row.names", all.x=TRUE)
rownames(annotated_results_up) <- annotated_results_up$Row.names
annotated_results_up$Row.names <- NULL
detected_up_NOISeq <- rownames(annotated_results_up)


#from differentially expressed genes, actually downregulated according to metadata
meta_up <- rownames(annotated_results[annotated_results$upregulation == 1,])

common_up <- intersect(detected_up_NOISeq, meta_up)

# Step 10: Summarize outliers
outliers_up <- setdiff(detected_up_NOISeq, meta_up)

# Step 11: Accuracy and Precision Matrix
true_positives <- length(common_up) + length(common_down)
false_positives <- length(setdiff(detected_up_NOISeq, meta_up)) + length(setdiff(detected_down_NOISeq, meta_down))
true_negatives <- nrow(meta_data) - (length(meta_up) + length(meta_down)) - false_positives
false_negatives <- length(setdiff(meta_up, detected_up_NOISeq)) + length(setdiff(meta_down, detected_down_NOISeq))
accuracy <- (true_positives + true_negatives) / nrow(meta_data)
precision <- true_positives / (true_positives + false_positives)
recall <- true_positives / (true_positives + false_negatives)
f1_score <- 2 * ((precision * recall) / (precision + recall))
FDR <- false_positives/(true_positives+false_positives) ## Precision = 1- FDR 

# Step 12: Output Metrics to CSV for further analysis outside R
metrics_df_5 <- data.frame(
  Threshold = q_value,
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

## Executing the same for q=0.4
q_value <- 0.4
DEgenes <- degenes(results, q = q_value, M = NULL)  

annotated_results <- merge(as.data.frame(DEgenes), meta_data, by="row.names", all.x=TRUE)
annotated_results
rownames(annotated_results) <- annotated_results$Row.names
annotated_results$Row.names <- NULL

# Downregulated in condition 2 compared to condition 1
detected_down = degenes(results, q = q_value, M = "up")
detected_down

annotated_results_down <- merge(as.data.frame(detected_down), meta_data, by="row.names", all.x=TRUE)
annotated_results_down
rownames(annotated_results_down) <- annotated_results_down$Row.names
annotated_results_down$Row.names <- NULL
detected_down_NOISeq <- rownames(annotated_results_down)
detected_down_NOISeq

#from differentially expressed genes, actually downregulated according to metadata
meta_down <- rownames(annotated_results[annotated_results$downregulation == 1,])

common_down <- intersect(detected_down_NOISeq, meta_down)

#detected_up represent the up regulated genes in condition 2 compared to condition 1.
detected_up = degenes(results, q = q_value, M = "down")

annotated_results_up <- merge(as.data.frame(detected_up), meta_data, by="row.names", all.x=TRUE)
rownames(annotated_results_up) <- annotated_results_up$Row.names
annotated_results_up$Row.names <- NULL
detected_up_NOISeq <- rownames(annotated_results_up)

#from differentially expressed genes, actually downregulated according to metadata
meta_up <- rownames(annotated_results[annotated_results$upregulation == 1,])
common_up <- intersect(detected_up_NOISeq, meta_up)
outliers_up <- setdiff(detected_up_NOISeq, meta_up)

true_positives <- length(common_up) + length(common_down)
false_positives <- length(setdiff(detected_up_NOISeq, meta_up)) + length(setdiff(detected_down_NOISeq, meta_down))
true_negatives <- nrow(meta_data) - (length(meta_up) + length(meta_down)) - false_positives
false_negatives <- length(setdiff(meta_up, detected_up_NOISeq)) + length(setdiff(meta_down, detected_down_NOISeq))
accuracy <- (true_positives + true_negatives) / nrow(meta_data)
precision <- true_positives / (true_positives + false_positives)
recall <- true_positives / (true_positives + false_negatives)
f1_score <- 2 * ((precision * recall) / (precision + recall))
FDR <- false_positives/(true_positives+false_positives) ## Precision = 1- FDR 

# Step 12: Output Metrics to CSV for further analysis outside R
metrics_df_6 <- data.frame(
  Threshold = q_value,
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

## Executing the same for q=0.3
q_value <- 0.3
DEgenes <- degenes(results, q = q_value, M = NULL)  

annotated_results <- merge(as.data.frame(DEgenes), meta_data, by="row.names", all.x=TRUE)
annotated_results
rownames(annotated_results) <- annotated_results$Row.names
annotated_results$Row.names <- NULL

# Downregulated in condition 2 compared to condition 1
detected_down = degenes(results, q = q_value, M = "up")
detected_down

annotated_results_down <- merge(as.data.frame(detected_down), meta_data, by="row.names", all.x=TRUE)
annotated_results_down
rownames(annotated_results_down) <- annotated_results_down$Row.names
annotated_results_down$Row.names <- NULL
detected_down_NOISeq <- rownames(annotated_results_down)
detected_down_NOISeq

#from differentially expressed genes, actually downregulated according to metadata
meta_down <- rownames(annotated_results[annotated_results$downregulation == 1,])

common_down <- intersect(detected_down_NOISeq, meta_down)

#detected_up represent the up regulated genes in condition 2 compared to condition 1.
detected_up = degenes(results, q = q_value, M = "down")

annotated_results_up <- merge(as.data.frame(detected_up), meta_data, by="row.names", all.x=TRUE)
rownames(annotated_results_up) <- annotated_results_up$Row.names
annotated_results_up$Row.names <- NULL
detected_up_NOISeq <- rownames(annotated_results_up)


#from differentially expressed genes, actually downregulated according to metadata
meta_up <- rownames(annotated_results[annotated_results$upregulation == 1,])

common_up <- intersect(detected_up_NOISeq, meta_up)

# Step 10: Summarize outliers
outliers_up <- setdiff(detected_up_NOISeq, meta_up)

# Step 11: Accuracy and Precision Matrix
true_positives <- length(common_up) + length(common_down)
false_positives <- length(setdiff(detected_up_NOISeq, meta_up)) + length(setdiff(detected_down_NOISeq, meta_down))
true_negatives <- nrow(meta_data) - (length(meta_up) + length(meta_down)) - false_positives
false_negatives <- length(setdiff(meta_up, detected_up_NOISeq)) + length(setdiff(meta_down, detected_down_NOISeq))
accuracy <- (true_positives + true_negatives) / nrow(meta_data)
precision <- true_positives / (true_positives + false_positives)
recall <- true_positives / (true_positives + false_negatives)
f1_score <- 2 * ((precision * recall) / (precision + recall))
FDR <- false_positives/(true_positives+false_positives) ## Precision = 1- FDR 

# Step 12: Output Metrics to CSV for further analysis outside R
metrics_df_7 <- data.frame(
  Threshold = q_value,
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
combined_data <- rbind(metrics_df_1, metrics_df_2,metrics_df_3,metrics_df_4,metrics_df_5,metrics_df_6,metrics_df_7)
combined_data
# Write the combined data frame back to the CSV file with append = TRUE
write.csv(combined_data, "Working Directory/Output/Threshold_noiseq_3_1000_0.csv", row.names = FALSE)

