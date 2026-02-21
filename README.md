*Yeast Genome Analysis Pipeline*
Scalable R pipeline for parsing multi-chromosome yeast genomic sequences, computing nucleotide composition, quantifying GC content, simulating gene expression levels, and performing statistical analysis with publication-quality visualizations.

Built using Biostrings and ggplot2.


**Project Overview**
Genomic sequence composition varies across chromosomes and can influence structural and functional properties such as transcriptional activity.

This project:
- Parses yeast chromosome sequence files
- Computes nucleotide frequencies (A, T, G, C)
- Calculates GC content per chromosome
- Simulates gene expression levels (for modeling demonstration)
- Performs statistical testing across chromosomes
- Generates error bar plots and regression analysis visualizations

The pipeline is modular and scalable for larger genomic datasets.

**Tools & Technologies**
- R
- Biostrings (Bioconductor) — sequence parsing & nucleotide frequency analysis
- ggplot2 — statistical visualization
- Base R statistical testing (one-way ANOVA equivalent, linear modeling)

**Pipeline Workflow**
1.Sequence Parsing
- Reads multiple chromosome .txt sequence files
- - Cleans and concatenates sequence lines
- Converts sequences to uppercase for consistency
- Uses alphabetFrequency() from Biostrings to compute base counts

Outputs per chromosome:

A count
T count
G count
C count
GC content (%)

2. GC Content Calculation
- GC content is computed for each chromosome and stored in a structured data frame.

3. Simulated Gene Expression
For modeling demonstration purposes:
 - Gene expression values are simulated using a normal distribution (mean = 10, sd = 2)
- Reproducibility ensured with set.seed(123)

4. Statistical Testing
To evaluate whether chromosome-level differences exist:
- One-way Welch’s ANOVA (oneway.test)
- Applied to:
      - GC_Content
      - Expression_Level

This tests whether mean differences between chromosomes are statistically significant.

