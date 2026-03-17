library(DESeq2)
library(topGO)
library(edgeR)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(stringr)

setwd('~/Documents/StA/popgen/scripts/Infestation_gene_expression/')

gaf<-read.table("~/Documents/StA/RNA/makeshift_gaf.tsv",h=F,sep='\t')
#create gene ontology overrep. functions
run_topGO<-function(res.df){
  gene_list<-res.df$padj
  names(gene_list)<-rownames(res.df)
  gene_list<-gene_list[!is.na(gene_list)]
  selection_function<-function(x) {
    return(x<0.05)#test geneset defined as Padj < 0.05; no logFC threshold
  }
  gene2GO<-list()
  for(i in 1:nrow(gaf)) {
    gene_symbol<-gaf[i,1] 
    go_terms<-gaf[i,2] 
    
    if(!is.na(go_terms) & !go_terms=="") {
      go_list<-trimws(unlist(strsplit(go_terms,";")))
      go_list<-go_list[go_list!=""]
      gene2GO[[gene_symbol]]<-go_list
    }
  }
  #test biological processes
  GOdata_BP <-  new("topGOdata",
                    ontology = "BP",
                    allGenes = gene_list,
                    geneSel = selection_function,
                    annot = annFUN.gene2GO,
                    gene2GO = gene2GO)
  #test molecular functions
  GOdata_MF <- new("topGOdata",
                   ontology = "MF",
                   allGenes = gene_list,
                   geneSel = selection_function,
                   annot = annFUN.gene2GO,
                   gene2GO = gene2GO)
  
  #run fisher tests
  test_BP <- runTest(GOdata_BP, algorithm = "weight01", statistic = "fisher")
  test_MF <- runTest(GOdata_MF, algorithm = "weight01", statistic = "fisher")
  results_BP <- GenTable(GOdata_BP, Fisher = test_BP, topNodes = 100)
  results_MF <- GenTable(GOdata_MF, Fisher = test_MF, topNodes = 100)
  results_BP$Category <- "Biological Process"
  results_MF$Category <- "Molecular Function"
  
  all_results <- rbind(results_BP, results_MF)
  all_results$Fisher<-as.numeric(all_results$Fisher)
  all_results<-all_results[all_results$Fisher<0.01,]#output only GO terms with P < 0.01
  return(all_results)
}

setwd('~/Documents/StA/RNA')
pheno<-read.csv('Sikkink_pheno.csv')
rownames(pheno)<-pheno$ID
cts<-read.csv('Sikkink_gene_count_matrix.csv',h=T,row.names="gene_id")
summary(colnames(cts) %in% rownames(pheno))
cts<-cts[,order(colnames(cts))]
pheno<-pheno[order(rownames(pheno)),]
summary(colnames(cts) == rownames(pheno))#check samples in same order

#DE analysis
dds <- DESeq(DESeqDataSetFromMatrix(countData = cts,
  colData = pheno,
  design= ~ Pop+Trt))#nb. two populations in analysis
resultsNames(dds)
summary(results(dds, name="Trt_D4_vs_C",alpha=0.05))#d4 vs c is the focus as it represents mid-stage infection
summary(results(dds, name="Trt_D7_vs_C",alpha=0.01))
res.d4<-data.frame(results(dds, name="Trt_D4_vs_C",alpha=0.05))
res.d4$gene<-rownames(res.d4)
res.d7<-data.frame(results(dds, name="Trt_D7_vs_C",alpha=0.05))
res.d7$gene<-rownames(res.d7)

#GO overrep
res.d4.go<-run_topGO(res.d4)
res.d4.go$fold<-res.d4.go$Significant/res.d4.go$Annotated
View(res.d4.go)


#plot volcano and GO overrep
g.volcano<-ggplot(res.d4,aes(x=log2FoldChange,y=-log10(pvalue)))+
  theme_bw()+theme(panel.grid=element_blank(),legend.position='none')+
  geom_point(size=0.75,aes(colour=padj<0.05))+
  geom_vline(xintercept=0,linetype='dashed')+
  xlim(c(-12,12))+
  scale_colour_manual(values=c("#aaaaaa","#C74955"))
g.go<-ggplot(res.d4.go[res.d4.go$Fisher<3.5e-5,],aes(y=Term,x=-log10(Fisher)))+
  theme_minimal()+theme(panel.grid=element_blank(),legend.position='none',
                        axis.title.y=element_blank(),axis.line.x=element_line())+
  geom_col(colour='black')+xlab("Fold enrichment")

g.volcano+g.go
ggsave('Infest_d4_vol_go.png',plot=g.volcano+labs(tag="F",title="Response to infestation")+
         g.go+theme(axis.text.y=element_text(size=8))+labs(tag="G")+
         plot_layout(widths=c(1.25,1)),height=2.6,width=7.5)


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
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] stringr_1.5.1               RColorBrewer_1.1-3          ggplot2_4.0.2              
# [4] pheatmap_1.0.12             edgeR_4.2.2                 limma_3.60.6               
# [7] topGO_2.56.0                SparseM_1.84-2              GO.db_3.19.1               
# [10] AnnotationDbi_1.66.0        graph_1.82.0                DESeq2_1.44.0              
# [13] SummarizedExperiment_1.34.0 Biobase_2.64.0              MatrixGenerics_1.16.0      
# [16] matrixStats_1.5.0           GenomicRanges_1.56.2        GenomeInfoDb_1.40.1        
# [19] IRanges_2.38.1              S4Vectors_0.42.1            BiocGenerics_0.50.0        
# 
# loaded via a namespace (and not attached):
#   [1] KEGGREST_1.44.1         gtable_0.3.6            lattice_0.22-6          vctrs_0.6.5            
# [5] tools_4.4.2             bitops_1.0-9            generics_0.1.3          parallel_4.4.2         
# [9] tibble_3.2.1            RSQLite_2.3.9           blob_1.2.4              pkgconfig_2.0.3        
# [13] Matrix_1.7-2            S7_0.2.1                lifecycle_1.0.4         GenomeInfoDbData_1.2.12
# [17] compiler_4.4.2          farver_2.1.2            Biostrings_2.72.1       statmod_1.5.0          
# [21] codetools_0.2-20        RCurl_1.98-1.16         pillar_1.10.1           crayon_1.5.3           
# [25] BiocParallel_1.38.0     DelayedArray_0.30.1     cachem_1.1.0            abind_1.4-8            
# [29] tidyselect_1.2.1        locfit_1.5-9.11         stringi_1.8.4           dplyr_1.1.4            
# [33] fastmap_1.2.0           grid_4.4.2              colorspace_2.1-1        cli_3.6.3              
# [37] SparseArray_1.4.8       magrittr_2.0.3          S4Arrays_1.4.1          withr_3.0.2            
# [41] scales_1.4.0            UCSC.utils_1.0.0        bit64_4.6.0-1           XVector_0.44.0         
# [45] httr_1.4.7              bit_4.5.0.1             png_0.1-8               memoise_2.0.1          
# [49] rlang_1.1.5             Rcpp_1.1.0              glue_1.8.0              DBI_1.2.3              
# [53] rstudioapi_0.17.1       jsonlite_1.8.9          R6_2.5.1                zlibbioc_1.50.0 