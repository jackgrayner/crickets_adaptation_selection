library(snpStats)
library(LDheatmap)
library(vcfR)

#define a function to convert the value into 0,1,2
convertToNumeric <- function(x){
  gdat <- matrix(NA,nrow = nrow(x), ncol = ncol(x))
  for (m in 1:nrow(x)){
    for (n in 1:ncol(x)){
      a <-as.numeric(unlist(strsplit(x[m,n], "|"))[1]) 
      
      b <- as.numeric(unlist(strsplit(x[m,n], "|"))[3])
      gdat[m,n] <- a+b
    }
  }
  rownames(gdat) <- rownames(x)
  colnames(gdat) <- colnames(x)
  return(gdat)
}

chr1=1
run_ld<-function(chr1){
  chr<-paste0('./ld_files/scaffold_',chr1)
  snp <- read.vcfR(paste0("LD_lostruct_mgeno0.01_",chr,".vcf.gz"))
  snpgt<-snp@gt
  snpgt_eur <- t(snpgt)
  #convert to snpMatrix - EUR
  snpgt_eur <- convertToNumeric(snpgt_eur)
  
  #load the snp_id_dist.csv, which contains the SNPs id and distance
  info <- read.table("LD_lostruct_mgeno0.01_200kb.map") 
  info<-info[info$V1==chr,]
  snpNames <- info$V2
  colnames(snpgt_eur) <- snpNames
  gdat_eur<-as(snpgt_eur,"SnpMatrix")
  
  mycols <- colorRampPalette(c("#C74955","#f79143","#555555" ))(30)
  
  ld.temp<-LDheatmap(gdat_eur,info$V4,add.map=FALSE,color = mycols,title = NULL,
                     #SNP.name = "scaffold_1:168831320:T:A",
                     add.key = TRUE,distances = "physical")
  ld<-ggdraw(ld.temp$LDheatmapGrob)
  return(ld)
}

LDheatmap(gdat_eur,info$V4,add.map=FALSE,color = mycols,title = NULL,
          add.key = TRUE,distances = "physical")


for (chr1 in c(1:14)){
  assign(paste0("LD",chr1),run_ld(chr1))
}


ggsave('ld_allchrs.png',plot=LD1+LD2+LD3+LD4+LD5+LD6+LD7+LD8+LD9+LD10+LD11+LD12+LD13+LD14+plot_layout(nrow=1),height=3,width=8)


# R version 4.4.2 (2024-10-31)
# Platform: aarch64-apple-darwin20
# Running under: macOS Sonoma 14.3
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
# LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# time zone: America/New_York
# tzcode source: internal
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] vcfR_1.15.0     LDheatmap_1.0-5 snpStats_1.54.0 Matrix_1.7-2    survival_3.8-3 
# 
# loaded via a namespace (and not attached):
#   [1] digest_0.6.37       permute_0.9-7       zlibbioc_1.50.0     mgcv_1.9-1          lattice_0.22-6     
# [6] magrittr_2.0.3      splines_4.4.2       parallel_4.4.2      BiocGenerics_0.50.0 viridisLite_0.4.2  
# [11] ape_5.8-1           pinfsc50_1.3.0      grid_4.4.2          vegan_2.6-10        compiler_4.4.2     
# [16] rstudioapi_0.17.1   tools_4.4.2         cluster_2.1.8       nlme_3.1-167        Rcpp_1.1.0         
# [21] MASS_7.3-64     