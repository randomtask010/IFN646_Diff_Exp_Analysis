run_loop_deseq2 <- function(SourceFileVariable, PValue) {
  
  library(DESeq2)
  library(ggplot2)
  
  # Load Dataset
  count_data <- read.table(paste0("RAW data/", SourceFileVariable, ".tsv"), header=TRUE, row.names=1)
  samplesize <- ncol(count_data) / 2
  sample_info <- data.frame(groups = factor(rep(1:2, each=samplesize)))
  
  # Create DESeqDataSet object
  dds <- DESeqDataSetFromMatrix(countData = count_data,
                                colData = sample_info,
                                design = ~ groups)
  
  # Filter low count genes
  dds <- dds[ rowSums(counts(dds)) >= 10, ]
  
  # Normalization
  dds <- estimateSizeFactors(dds)
  normalized_counts <- counts(dds, normalized=TRUE)
  
  # Perform differential expression analysis
  dds <- DESeq(dds)
  
  # Extract results
  res <- results(dds)
  summary(res)
  res <- res[!is.na(res$padj < PValue), ]
  
  # Load Metadata
  meta_data <- read.table(paste0("RAW data/", SourceFileVariable, "_meta.tsv"), header=TRUE, row.names=1)
  
  # Merge DESeq2 results with metadata
  annotated_results <- merge(as.data.frame(res), meta_data, by="row.names", all.x=TRUE)
  rownames(annotated_results) <- annotated_results$Row.names
  annotated_results$Row.names <- NULL
  
  # Comparison using LogFC
  detected_up_DESeq2 <- rownames(annotated_results[annotated_results$log2FoldChange > 1 & annotated_results$padj < PValue,])
  meta_up <- rownames(annotated_results[annotated_results$upregulation == 1,])
  detected_down_DESeq2 <- rownames(annotated_results[annotated_results$log2FoldChange < -1 & annotated_results$padj < PValue,])
  meta_down <- rownames(annotated_results[annotated_results$downregulation == 1,])
  
  # Identify false positives
  outliers_up <- setdiff(detected_up_DESeq2, meta_up)
  outliers_down <- setdiff(detected_down_DESeq2, meta_down)
  false_positives_genes <- c(outliers_up, outliers_down)
  
  # Filter results for false positives
  res_fp <- as.data.frame(res[rownames(res) %in% false_positives_genes, ])
  
  
  # Volcano plot for false positives
  volcano_plot_fp <- ggplot(res_fp, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = padj < PValue & abs(log2FoldChange) > 1), alpha = 0.7) +
    geom_hline(yintercept = -log10(PValue), linetype = "dashed") +   # Update here
    labs(
      title = paste("Volcano Plot for False Positives in", SourceFileVariable),
      x = "Log2 Fold Change",
      y = "-Log10 Adjusted P-value"
    ) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
    theme_bw()
  
  
  print(volcano_plot_fp)
  
  # Save the plot
  output_file_volcano_fp <- paste0("Working Directory/Output/Images/Volcano_FP_DESeq2", SourceFileVariable, ".png")
  ggsave(output_file_volcano_fp, plot = volcano_plot_fp, width = 10, height = 6)
}
