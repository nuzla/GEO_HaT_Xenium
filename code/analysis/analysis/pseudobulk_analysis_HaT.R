
# Load libraries
library(Seurat)
library(DESeq2)
library(Matrix)
library(ggplot2)
library(ggrepel)
library(dplyr)

# Sample metadata table
sample_info <- data.frame(
    sample_id = c("18-339", "17-476", "18-031", "18-201", "17-432", "17-044", "18-070", "18-547"),
    hat_status = c("HaT", "HaT", "HaT", "HaT", "non-HaT", "non-HaT", "non-HaT", "non-HaT"),
    expr_path = c("path/to/18-339/cell_feature_matrix.h5", "path/to/17-476/cell_feature_matrix.h5", 
                  "path/to/18-031/cell_feature_matrix.h5", "path/to/18-201/cell_feature_matrix.h5", 
                  "path/to/17-432/cell_feature_matrix.h5", "path/to/17-044/cell_feature_matrix.h5", 
                  "path/to/18-070/cell_feature_matrix.h5", "path/to/18-547/cell_feature_matrix.h5"),
    meta_path = c("path/to/18-339/cells.csv.gz", "path/to/17-476/cells.csv.gz", 
                  "path/to/18-031/cells.csv.gz", "path/to/18-201/cells.csv.gz", 
                  "path/to/17-432/cells.csv.gz", "path/to/17-044/cells.csv.gz", 
                  "path/to/18-070/cells.csv.gz", "path/to/18-547/cells.csv.gz")
)

# Load function
load_sample <- function(expr_path, meta_path, sample_id) {
  seurat_obj <- LoadXenium(expr_path)
  metadata <- read.csv(meta_path)
  metadata$barcode <- make.unique(as.character(metadata$cell_id))
  rownames(metadata) <- metadata$barcode
  metadata$sample_id <- sample_id
  seurat_obj <- AddMetaData(seurat_obj, metadata)
  seurat_obj$sample_id <- sample_id
  return(seurat_obj)
}

# Load samples
seurat_list <- mapply(load_sample,
                      expr_path = sample_info$expr_path,
                      meta_path = sample_info$meta_path,
                      sample_id = sample_info$sample_id,
                      SIMPLIFY = FALSE)

# Merge
combined <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = sample_info$sample_id)

# Extract and combine expression matrices
expr_layers <- grep("counts.Gene Expression", Layers(combined[["RNA"]]), value = TRUE)
counts_list <- lapply(expr_layers, function(layer) {
    mat <- combined[["RNA"]][[layer]]
    if (is.data.frame(mat) && ncol(mat) > 0 && nrow(mat) > 0) {
        mat <- as.matrix(mat)
        rownames(mat) <- make.unique(rownames(mat))
        return(mat)
    } else {
        return(NULL)
    }
})
counts_list <- counts_list[!sapply(counts_list, is.null)]
counts <- Reduce("+", counts_list)

# Pseudobulk aggregation
meta <- combined@meta.data
meta$sample_id <- factor(meta$sample_id, levels = sample_info$sample_id)
pseudobulk_counts_list <- lapply(levels(meta$sample_id), function(sid) {
    cell_ids <- rownames(meta)[meta$sample_id == sid]
    rowSums(counts[, cell_ids, drop = FALSE])
})
names(pseudobulk_counts_list) <- levels(meta$sample_id)
pseudobulk_counts <- do.call(cbind, pseudobulk_counts_list)

# Create colData
sample_table <- data.frame(
    row.names = colnames(pseudobulk_counts),
    group = sample_info$hat_status
)

# DESeq2
dds <- DESeqDataSetFromMatrix(countData = pseudobulk_counts,
                              colData = sample_table,
                              design = ~ group)
dds <- DESeq(dds)
res <- results(dds, contrast = c("group", "HaT", "non-HaT"))

# Volcano plot
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df <- res_df[complete.cases(res_df), ]
res_df$significance <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, "significant", "not_significant")
res_df$highlight <- ifelse(res_df$gene %in% c("MRGPRX2", "TPSAB1", "LAMP1", "CPA3", "ENPP3", "CD63", "SIGLEC8", "BAFT", "ASCL2"), TRUE, FALSE)

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = significance), alpha = 0.5) +
  geom_text_repel(data = subset(res_df, highlight),
                  aes(label = gene),
                  size = 3, max.overlaps = Inf) +
  theme_minimal() +
  scale_color_manual(values = c("gray", "red")) +
  ggtitle("Volcano Plot - HaT vs non-HaT")

# Boxplot for TPSAB1
gene <- "TPSAB1"
expr_vals <- counts[gene, ]
box_df <- data.frame(
  Expression = expr_vals,
  Sample = names(expr_vals),
  Group = sample_info$hat_status[match(names(expr_vals), sample_info$sample_id)]
)

ggplot(box_df, aes(x = Group, y = Expression, fill = Group)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  labs(title = paste("Expression of", gene)) +
  theme_minimal()
