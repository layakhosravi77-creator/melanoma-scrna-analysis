
# GSE72056: Initial melanoma scRNA-seq analysis
# Dataset: GSE72056
# Purpose: Explore immune-related gene expression
# Note: Confirm expression scale before further normalization.

# 1. Load the downloaded dataset
file_path <- "data/GSE72056_melanoma_single_cell_revised_v2.txt.gz"

geo_data <- read.delim(
  gzfile(file_path),
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  quote = "",
  stringsAsFactors = FALSE
)

# 2. Extract cell annotations
cell_metadata <- data.frame(
  tumor = as.character(unlist(geo_data[1, -1], use.names = FALSE)),
  malignant = as.integer(unlist(geo_data[2, -1], use.names = FALSE)),
  non_malignant_type = as.integer(
    unlist(geo_data[3, -1], use.names = FALSE)
  ),
  row.names = colnames(geo_data)[-1],
  check.names = FALSE
)

# 3. Extract gene expression matrix
expression_df <- geo_data[-(1:3), ]
gene_names <- as.character(expression_df[[1]])

expression_matrix <- as.matrix(expression_df[, -1])
storage.mode(expression_matrix) <- "numeric"
rownames(expression_matrix) <- make.unique(gene_names)

# Release the intermediate data frame
rm(expression_df, geo_data)

# 4. Check dimensions and expression range
print(dim(expression_matrix))
print(dim(cell_metadata))
print(range(expression_matrix, na.rm = TRUE))

# 5. Summarize cell annotations
print(table(cell_metadata$malignant, useNA = "ifany"))
print(table(cell_metadata$non_malignant_type, useNA = "ifany"))

# 6. Select T cells and malignant cells
t_cells <- cell_metadata$malignant == 1 &
  cell_metadata$non_malignant_type == 1

malignant_cells <- cell_metadata$malignant == 2

print(sum(t_cells))
print(sum(malignant_cells))

# 7. Select genes related to immune regulation
genes_of_interest <- c(
  "PDCD1", "CD274", "PDCD1LG2",
  "LAG3", "HAVCR2", "TIGIT",
  "TOX", "CTLA4", "GZMB", "IFNG"
)

genes_found <- intersect(
  genes_of_interest,
  rownames(expression_matrix)
)

print(genes_found)

# 8. Compare mean expression and detected cells
gene_comparison <- data.frame(
  gene = genes_found,
  T_cell_mean = rowMeans(
    expression_matrix[genes_found, t_cells, drop = FALSE]
  ),
  malignant_mean = rowMeans(
    expression_matrix[genes_found, malignant_cells, drop = FALSE]
  ),
  T_cells_detected = rowSums(
    expression_matrix[genes_found, t_cells, drop = FALSE] > 0
  ),
  malignant_cells_detected = rowSums(
    expression_matrix[genes_found, malignant_cells, drop = FALSE] > 0
  )
)

print(gene_comparison)

# 9. Plot gene expression in T cells
par(mar = c(8, 4, 3, 1))
barplot(
  gene_comparison$T_cell_mean,
  names.arg = gene_comparison$gene,
  las = 2,
  main = "Gene expression in T cells",
  ylab = "Mean expression",
  col = "skyblue",
  cex.names = 0.8
)

# 10. Plot gene expression in malignant cells
par(mar = c(8, 4, 3, 1))
barplot(
  gene_comparison$malignant_mean,
  names.arg = gene_comparison$gene,
  las = 2,
  main = "Gene expression in malignant cells",
  ylab = "Mean expression",
  col = "salmon",
  cex.names = 0.8
)
