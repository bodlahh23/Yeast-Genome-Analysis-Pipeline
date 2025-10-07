# Install necessary packages if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  BiocManager::install("Biostrings")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

# Load libraries
library(Biostrings)
library(ggplot2)

# Define the folder path where all yeast chromosome text files are located
folder_path <- "C:/Users/bodla/Downloads/course-code-repo-main/course-code-repo-main/bioinformatics_(BIOL130-230)/genomics_gene-parsing-processing/YeastGenes"

# Get list of all text files
file_list <- list.files(path = folder_path, pattern = "*.txt", full.names = TRUE)

# Function to read each text file and extract base counts and GC content
process_file <- function(file_path) {
  seq_lines <- readLines(file_path)
  # Remove possible empty lines
  seq_lines <- seq_lines[seq_lines != ""]
  sequence <- toupper(paste(seq_lines, collapse = ""))
  
  # Count nucleotides
  counts <- alphabetFrequency(DNAString(sequence), baseOnly = TRUE)
  total_bases <- sum(counts[c("A", "T", "G", "C")])
  gc_content <- 100 * (counts["G"] + counts["C"]) / total_bases
  
  # Extract Chromosome info from filename
  file_name <- basename(file_path)
  chrom_name <- tools::file_path_sans_ext(file_name)
  
  data.frame(
    Chromosome = chrom_name,
    A = counts["A"],
    T = counts["T"],
    G = counts["G"],
    C = counts["C"],
    GC_Content = gc_content
  )
}

# Apply function to all files and create one big data frame
nuc_data <- do.call(rbind, lapply(file_list, process_file))

# Simulate Expression Levels (replace with real data if you have)
set.seed(123)  # For reproducibility
nuc_data$Expression_Level <- rnorm(nrow(nuc_data), mean = 10, sd = 2)

# Make Chromosome a factor
nuc_data$Chromosome <- as.factor(nuc_data$Chromosome)

# Save a CSV for your report if needed
write.csv(nuc_data, "yeast_chromosome_summary.csv", row.names = FALSE)

# ---- Statistical Testing ----

# Function to perform oneway.test
test_metric <- function(metric_name) {
  formula <- as.formula(paste(metric_name, "~ Chromosome"))
  result <- oneway.test(formula, data = nuc_data)
  cat("Oneway Test for", metric_name, ":\n")
  print(result)
  cat("\n")
}

# Perform oneway tests
test_metric("GC_Content")
test_metric("Expression_Level")

# ---- Error Bar Plots ----

plot_error_bars <- function(metric_name) {
  summary_data <- aggregate(nuc_data[[metric_name]],
                            by = list(Chromosome = nuc_data$Chromosome),
                            FUN = function(x) c(mean = mean(x), sd = sd(x)))
  summary_data <- do.call(data.frame, summary_data)
  colnames(summary_data)[2:3] <- c("Mean", "SD")
  
  p <- ggplot(summary_data, aes(x = Chromosome, y = Mean)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.3) +
    ylab(metric_name) +
    ggtitle(paste("Error Bar Plot for", metric_name)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p)
}

# Generate error bar plots
plot_error_bars("GC_Content")
plot_error_bars("Expression_Level")

# ---- Scatterplot with Linear Regression ----

scatterplot_with_regression <- function(x_metric, y_metric) {
  p <- ggplot(nuc_data, aes_string(x = x_metric, y = y_metric)) +
    geom_point(color = "darkred", alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, color = "black") +
    ggtitle(paste(x_metric, "vs", y_metric)) +
    theme_minimal()
  
  print(p)
  
  # Linear model summary
  model <- lm(as.formula(paste(y_metric, "~", x_metric)), data = nuc_data)
  cat("\nLinear Model Summary for", x_metric, "vs", y_metric, ":\n")
  print(summary(model))
}

# Example scatterplot
scatterplot_with_regression("GC_Content", "Expression_Level")

