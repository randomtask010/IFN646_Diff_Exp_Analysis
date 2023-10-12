run_loop_edgeR_threshold <-function(SourceFileVariable, PValue) {
  
  
  # Step 1: Call Library (edgeR)
  library(edgeR)
  
  # Step 2: Load Dataset
  count_data <- read.table(paste0("RAW data/", SourceFileVariable, ".tsv"), header=TRUE, row.names=1)
  samplesize <- ncol(count_data) /2
  groups <- factor(rep(1:2, each=samplesize))
  
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
  decided <- decideTestsDGE(et, adjust.method="BH", p.value=PValue)
  summary(decided)
  
  # Step 10: Print top genes
  topTags(et)
  
  # Step 11: Load Metadata
  meta_data <- read.table(paste0("RAW data/", SourceFileVariable,"_meta.tsv"), header=TRUE, row.names=1)
  edgeR_results <- topTags(et, n=Inf)
  annotated_results <- merge(edgeR_results, meta_data, by="row.names", all.x=TRUE)
  rownames(annotated_results) <- annotated_results$Row.names
  annotated_results$Row.names <- NULL
  
  # Step 12: Comparison using LogFC
  # Upregulated
  detected_up_edgeR <- rownames(annotated_results[annotated_results$logFC > 0 & annotated_results$FDR < PValue,])
  meta_up <- rownames(annotated_results[annotated_results$upregulation == 1,])
  common_up <- intersect(detected_up_edgeR, meta_up)
  # Downregulated
  detected_down_edgeR <- rownames(annotated_results[annotated_results$logFC < 0 & annotated_results$FDR < PValue,])
  meta_down <- rownames(annotated_results[annotated_results$downregulation == 1,])
  common_down <- intersect(detected_down_edgeR, meta_down)
  
  # Step 13: Summarize outliers
  outliers_up <- setdiff(detected_up_edgeR, meta_up)
  write.csv(outliers_up, paste0("Working Directory/Output/",Tool,"_" , SourceFileVariable,"_outliers_upregulated_",  "PValue_",PValue,".csv"), row.names = FALSE)
  outliers_down <- setdiff(detected_down_edgeR, meta_down)
  write.csv(outliers_down, paste0("Working Directory/Output/", Tool,"_", SourceFileVariable, "_outliers_downregulated_", "PValue_", PValue,".csv"), row.names = FALSE)
  
  # Step 14: Accuracy and Precision Matrix
  true_positives <- length(common_up) + length(common_down)
  false_positives <- length(setdiff(detected_up_edgeR, meta_up)) + length(setdiff(detected_down_edgeR, meta_down))
  true_negatives <- nrow(meta_data) - (length(meta_up) + length(meta_down)) - false_positives
  false_negatives <- length(setdiff(meta_up, detected_up_edgeR)) + length(setdiff(meta_down, detected_down_edgeR))
  accuracy <- (true_positives + true_negatives) / nrow(meta_data)
  precision <- true_positives / (true_positives + false_positives)
  recall <- true_positives / (true_positives + false_negatives)
  f1_score <- 2 * ((precision * recall) / (precision + recall))
  fdr <- false_positives / (false_positives + true_positives)
  
  # Step 15: Output Metrics to CSV for further analysis outside R
  metrics_df <- data.frame(
    Threshold = PValue,
    True_Positives = true_positives,
    False_Positives = false_positives,
    True_Negatives = true_negatives,
    False_Negatives = false_negatives,
    Accuracy = accuracy,
    Precision = precision,
    Recall = recall,
    F1_Score = f1_score,
    FDR = fdr,
    Experiment = SourceFileVariable
  )

  # Initialize an empty data frame to hold all the results
  all_results <- data.frame()
  
  # Read the existing results, if any
  if (file.exists("Working Directory/Output/Threshold_edger.csv")) {
    all_results <- read.csv("Working Directory/Output/Threshold_edger.csv", header = TRUE)
  }
  
  # Append the current results to the all_results data frame
  all_results <- rbind(all_results, metrics_df)
  
  # Write the combined data frame to the CSV file
  write.csv(all_results, "Working Directory/Output/Threshold_edger.csv", row.names = FALSE)
  
  
}