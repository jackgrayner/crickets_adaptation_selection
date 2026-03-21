#!/bin/bash
#SBATCH --job-name=bcf_merge_concat  
#SBATCH --cpus-per-task=16   
#SBATCH --mem=40gb                
#SBATCH --partition=long
#SBATCH --output=bcf_merge_concat.log
pwd; hostname; date

source activate bcftools
bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' ./Variants_all_chr.vcf.gz -O z -o Variants_all_chr_anno.vcf.gz
bcftools index Variants_all_chr_anno.vcf.gz

source activate vcftools
vcftools --gzvcf Variants_all_chr_anno.vcf.gz  --maf 0.1 --min-alleles 2 --max-alleles 2 --max-missing 0.1 --minDP 5 --minGQ 20 --recode --recode-INFO-all --stdout | bgzip -c > Variants_all_chr_anno_filt.vcf.gz

source activate bcftools
bcftools index Variants_all_chr_anno_filt.vcf.gz
#nb allpops_alloutgroup_rad_vars.vcf.gz is the full Hawaiian vcf subset for variants covered by the rad-seq data
bcftools merge Variants_all_chr_anno_filt.vcf.gz allpops_alloutgroup_rad_vars.vcf.gz -Oz -o RAD_WGS_merged.vcf.gz
bcftools index RAD_WGS_merged.vcf.gz

source activate vcftools
vcftools --gzvcf RAD_WGS_merged.vcf.gz  --maf 0.05 --min-alleles 2 --max-alleles 2 --max-missing 0.5 --minDP 5 --minGQ 20 --recode --recode-INFO-all --stdout | bgzip -c > RAD_WGS_merged_filt.vcf.gz

cd ./pca
source activate plink2
plink2 --allow-extra-chr --vcf ../RAD_WGS_merged_filt.vcf.gz \
 --double-id --geno 0.7 --max-alleles 2 --make-bed --maf 0.05 --snps-only \
 --out RAD_WGS_merged

for chr1 in {1..14}
do
	chr=scaffold_${chr1}
	plink2 --allow-extra-chr --remove exclude.txt --bfile RAD_WGS_merged \
		 --chr ${chr} --make-bed --mind 0.95 --snps-only --pca \
		 --out merged_rad_${chr}
done
