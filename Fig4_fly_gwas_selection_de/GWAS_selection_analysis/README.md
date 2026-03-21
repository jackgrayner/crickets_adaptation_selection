# GWAS and Selection analyses

## association_test.sh
Takes a bed/fam file as input, calculates a relatedness matrix (GEMMA), and runs a linear mixed model (GEMMA) to test for genoetype-phenotype associations, accounting for relatedness and population stratification.

## genome_selection_analysis.R
Reads in a ped/map file and a csv file with the predictor variable (fly attack rate, i.e., selection), as well as a K value (num. of genetic groups in data) passed when submitting the job. The ped file is converted to LFMM format and missing genotypes are imputed with LEA's impute() function before running a latent factor mixed model (LFMM) to test for selection-associated variants.

## plotGWAS_selection_results.R
This script:
1. Reads in the results files from the above analyses, and identifies overlapping and nearby genes.
2. Plots genetic associations for Fw on Chr1 and Cw on Chr2, based on prior knowledge of their genetic architecture from mapping studies.
3. Plots the strength of association between fly attack rates and genotype across the whole genome (see panels C - E in plotGWAS_selection_example.png), and tests for overrepresentation of gene ontology categories (topGO). 

Raw results files are too large to upload so I have uploaded preview files including the first 10,000 lines.

## Results files

The top 200 variants (ranked by P-value) for each of the GWAS are included: 
- Cw_top200.csv: Cw phenotype, all samples
- Fw_Kauai_top200.csv: Fw phenotype, Kauai samples only
- Fw_Oahu_top200.csv: Fw phenotype, Oahu samples only

Columns (except for last two) are from the standard GEMMA output:

- chr: Chromosome number
- rs: SNP ID 
- ps: Physical position in base pairs
- n_miss: Number of missing individuals for this SNP
- af: Minor allele frequency (MAF)
- logl_H1: log likelihood under the alternative hypothesis 
- l_mle: log-likelihood at the maximum likelihood estimate
- p_lrt: p-value for the Likelihood Ratio Test 
- Island (if present): island included in test
- Padj: BH-adjusted p_lrt
- sig: significant at Padj < 0.05 (N or Y)
