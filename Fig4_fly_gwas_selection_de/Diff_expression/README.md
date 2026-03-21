## Differential expression analyses/results for infestation data

Data from:
K. L. Sikkink, N. W. Bailey, M. Zuk, S. L. Balenger, Immunogenetic and tolerance strategies against a novel parasitoid of wild field crickets. Ecology and Evolution 10, 13312–13326 (2020).

### Files

- Sikkink_gene_count_matrix.csv - gene abundances per sample
- Sikkink_pheno.csv - phenotype information per sample
- DESeq2_infestation_results.csv 
    - baseMean: The average of normalized counts across all samples, accounting for library size.
    - log2FoldChange (LFC): The effect size estimate, showing how much a gene’s expression differs between groups (e.g., Treatment vs. Control) on a scale. Positive values indicate up-regulation in the numerator condition (e.g., treated), while negative values indicate down-regulation.
    - lfcSE: The standard error of the log2 fold change estimate.
    - stat: The Wald statistic (for Wald test) or the difference in deviance (for Likelihood Ratio Test)
    - pvalue: The Wald test -value for the significance of the log2 fold change.
    - padj: The P-value adjusted for multiple testing using the Benjamini-Hochberg method to control the False Discovery Rate (FDR).
- stringtie_annotation_positions.csv - gene locations from genome annotation
