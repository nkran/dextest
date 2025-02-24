
# Load necessary libraries
library(DESeq2)
library(dplyr)
library(tximport)
library(tibble)
library(ggplot2)
library(reshape2)
library(clusterProfiler)

# Define file paths
samples <- read.table("samples.txt", header = TRUE)
files <- file.path("kallisto_results", samples$sample_id, "abundance.tsv")
names(files) <- samples$run

txi <- tximport(files, type = "kallisto", txOut = TRUE)

dds <- DESeqDataSetFromTximport(
  txi, colData = samples, design = ~ batch_time
)

dds <- DESeq(dds)

for (i in 1:8) {
  contrast_name <- paste0("B", i, "_vs_A", i)
  res_ba <- results(
    dds, contrast = c("batch_time", paste0("B", i), paste0("A", i)), 
    alpha = 0.05, lfcThreshold = 1, altHypothesis = "greaterAbs"
  )
  res_ape <- lfcShrink(
    dds = dds, contrast = c("batch_time", paste0("B", i), paste0("A", i)), 
    res = res_ba, type = "normal"
  )
  df_ba <- as.data.frame(res_ba) %>% rownames_to_column(var = "gene_id")
  
  # Use ggplot2 to plot volcano plot, plot only statistically significant genes
  p <- ggplot(df_ba, aes(x = log2FoldChange, y = -log10(pvalue))) +
    geom_point(aes(color = ifelse(padj < 0.05, "red", "black"))) +
    theme_minimal() +
    theme(legend.position = "none") +
    xlab("log2 fold change") +
    ylab("-log10 p-value") +
    ggtitle(paste("Volcano plot of differentially expressed genes -", contrast_name)) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "red") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red")
  
  # Save plot
  ggsave(
    filename = paste0("results_dex/volcano_", contrast_name, ".png"), plot = p
  )
  
  # Save results
  write.csv(df_ba, paste0("results_dex/res_", contrast_name, ".csv"))
}
  write.csv(df.BA, paste0("results_dex/res_", contrast_name, ".csv"))
}




# vsd <- vst(dds, nsub=5000)
# # plot PCA
# # plotPCA(vsd, intgroup = c("timepoint", "condition"))
# plotPCA(vsd, intgroup = c("batch_time"))






# df.BA.sig <- df.BA %>% filter(padj<0.05) %>% arrange(log2FoldChange)


# norm.counts <- as.data.frame(counts(dds, normalized = TRUE))
# colnames(norm.counts) <- samples$sample_id
# norm.counts$gene_id <- rownames(norm.counts)
# norm.counts.sig <- norm.counts[df.BA.sig$gene_id,]

# norm.counts.sig.melt <- melt(norm.counts.sig, id.vars = c("gene_id"), variable.name = "sample_id", value.name = "normcounts")
# norm.counts.sig.melt <- merge(norm.counts.sig.melt, samples, by = "sample_id")


