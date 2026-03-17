library(lfmm)
library(LEA)

kval <- commandArgs(trailingOnly=TRUE)[1]#read kval from command line

#convert ped file to LFMM format
#input file - ped file pruned using plink's --indep-pairwise 50kb 1 0.1
ped2lfmm("allpops_nodupes_filtered_pruned.ped", 
         output.file = "allpops_nodupes_filtered_pruned.lfmm", force = TRUE)
Y<-read.lfmm('allpops_nodupes_filtered_pruned.lfmm')
pc <- prcomp(Y)
plot(pc$sdev[1:20]^2, xlab = 'PC', ylab = "Variance explained")
points(6,pc$sdev[6]^2, type = "h", lwd = 3, col = "blue")

project.snmf = snmf("allpops_nodupes_filtered_pruned.lfmm", K = kval, 
                                        entropy = TRUE, repetitions = 10,
                                        project = "new")

#impute missing genotypes
best = which.min(cross.entropy(project.snmf, K = kval))
impute(project.snmf, "allpops_nodupes_filtered_pruned.lfmm", method = 'mode', K = kval, run = best)

#png(paste0(chr,'_pca_scores.png'))
#plot(pc$sdev[1:20]^2, xlab = 'PC', ylab = "Variance explained")
#points(6,pc$sdev[6]^2, type = "h", lwd = 3, col = "blue")
#dev.off()

#read in imputed data
Y<-read.lfmm('allpops_nodupes_filtered_pruned.lfmm_imputed.lfmm')
#read in selection values
s=read.csv("selection.csv")
s1=s$selection

#run model
mod.lfmm <- lfmm_ridge(Y = Y, 
                       X = s1 , 
                       K = kval)
#test model
pv <- lfmm_test(Y = Y, 
                X = s1 , 
                lfmm = mod.lfmm, 
                calibrate = "gif")

pvalues <- pv$calibrated.pvalue 

# qqplot(rexp(length(pvalues), rate = log(10)),
#        -log10(pvalues), xlab = "Expected quantile",
#        pch = 19, cex = .4)
# abline(0,1)

# plot(-log10(pvalues), 
#      pch = 19, 
#      cex = .2, 
#      xlab = "SNP", ylab = "-Log P",
#      col = "grey")

#read in variant information and merge for final results file
map.file<-read.table("allpops_nodupes_filtered_pruned.map",h=F)
colnames(map.file)<-c("chr","variant","V3","ps")
map.file$P<-pvalues
map.file<-map.file[!is.na(map.file$P),]
map.file$Padj<-p.adjust(map.file$P)

write.csv(map.file,"LFMM_results.csv",row.names=FALSE,quote=FALSE)


# R version 4.4.2 (2024-10-31)
# Platform: x86_64-conda-linux-gnu
# Running under: Debian GNU/Linux 12 (bookworm)
# 
# Matrix products: default
# BLAS/LAPACK: /mnt/apps/users/jrayner/conda/lib/libopenblasp-r0.3.28.so;  LAPACK version 3.12.0
# 
# locale:
#   [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C              
# [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
# [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
# [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: Europe/London
# tzcode source: system (glibc)
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] LEA_3.18.0 lfmm_1.1  
# 
# loaded via a namespace (and not attached):
#   [1] compiler_4.4.2   Rcpp_1.0.14      codetools_0.2-20 iterators_1.0.14
# [5] foreach_1.5.2   