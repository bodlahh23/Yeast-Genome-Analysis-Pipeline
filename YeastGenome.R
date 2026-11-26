

# Install packages if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("Biostrings", quietly = TRUE)) {
  BiocManager::install("Biostrings")
}

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

library(Biostrings)
library(ggplot2)

folder_path <- "C:/Users/bodla/Downloads/course-code-repo-main/course-code-repo-main/bioinformatics_(BIOL130-230)/genomics_gene-parsing-processing/YeastGenes"

file_list <- list.files(
  path = folder_path,
  pattern = "\\.txt$",
  full.names = TRUE
)

process_file <- function(file_path) {
  
  
  seq_lines <- suppressWarnings(readLines(file_path))
  seq_lines <- seq_lines[nchar(seq_lines) > 0]
  seq_lines <- seq_lines[!grepl("^>", seq_lines)]
  sequence <- toupper(paste(seq_lines, collapse = ""))
  sequence <- gsub("[^ATGCN]", "", sequence)
  dna <- DNAString(sequence)
  counts <- alphabetFrequency(dna, baseOnly = TRUE)
  total_bases <- sum(counts[c("A", "T", "G", "C")])
  gc_content <- ((counts["G"] + counts["C"]) / total_bases) * 100
  gene_name <- tools::file_path_sans_ext(basename(file_path))
  
  data.frame(
    Gene = gene_name,
    Length = total_bases,
    A = counts["A"],
    T = counts["T"],
    G = counts["G"],
    C = counts["C"],
    GC_Content = gc_content
  )
}


nuc_data <- do.call(
  rbind,
  lapply(file_list, process_file)
)


head(nuc_data)

cat("\nNumber of genes processed:", nrow(nuc_data), "\n")


set.seed(123)

nuc_data$Expression_Level <- rnorm(
  nrow(nuc_data),
  mean = 10,
  sd = 2
)


write.csv(
  nuc_data,
  "yeast_gene_summary.csv",
  row.names = FALSE
)


cat("\nSummary Statistics\n")
print(summary(nuc_data))


ggplot(nuc_data, aes(x = GC_Content)) +
  geom_histogram(
    bins = 20,
    fill = "steelblue",
    color = "black"
  ) +
  theme_minimal() +
  labs(
    title = "Distribution of GC Content",
    x = "GC Content (%)",
    y = "Count"
  )

ggplot(nuc_data, aes(x = Length)) +
  geom_histogram(
    bins = 20,
    fill = "forestgreen",
    color = "black"
  ) +
  theme_minimal() +
  labs(
    title = "Distribution of Gene Lengths",
    x = "Length (bp)",
    y = "Count"
  )

ggplot(
  nuc_data,
  aes(x = GC_Content, y = Expression_Level)
) +
  geom_point(
    color = "darkred",
    alpha = 0.7
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "black"
  ) +
  theme_minimal() +
  labs(
    title = "GC Content vs Expression Level",
    x = "GC Content (%)",
    y = "Expression Level"
  )


model <- lm(
  Expression_Level ~ GC_Content,
  data = nuc_data
)

cat("\n\nLinear Model Results:\n")
print(summary(model))


cat("\n\nCorrelation Test:\n")
print(cor.test(
  nuc_data$GC_Content,
  nuc_data$Expression_Level
))


top_gc <- nuc_data[
  order(-nuc_data$GC_Content),
][1:10, ]

cat("\nTop 10 Genes by GC Content:\n")
print(top_gc)


write.csv(
  top_gc,
  "top_gc_genes.csv",
  row.names = FALSE
)
