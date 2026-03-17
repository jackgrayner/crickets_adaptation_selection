library(ggtree)
library(ggplot2)
library(dplyr)
library(ape)
library(ggExtra)
library(patchwork)
library(pcadapt)
library(ggplotify)
library(tidyr)
library(tidyverse)
library(ggbeeswarm)
library(plyr)
library(ggplot2)
library(dplyr)
library(patchwork)
library(pheatmap)

science_theme <- theme_bw(base_size = 6) +
  theme(
    panel.grid=element_blank(),
    text = element_text(family = "Arial"),
    plot.title = element_text(size = 10),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7.5),
    strip.text = element_text(size = 8),
    plot.tag = element_text(size = 9, face = "bold"),
    plot.background=element_rect(fill='white',colour='white'),
    legend.position='none'
  )

#######
# PCA
#######

info <- read.table("./pca_pops/allHawaii_nodupes_pruned.fam",h=F)
pops<-info$V1
filename <- read.pcadapt("./pca_pops/allpops_outgroup_nodupes_chr5_9_14_pruned.bed", type = "bed")
auto <- pcadapt(input = filename, K = 9,LD.clumping = list(size = 500, thr = 0.1))
plot(auto, option = "screeplot")#k=9 looks reasonable
pca.df<-data.frame(auto$scores)
pca.df<-pca.df[c(1:300),]
pca.df$population<-info$V1
pca.df$ID<-info$V2
pca.df$ID <- sub(".*-", "", pca.df$ID)
pca.df$island <- sub("_.*", "", pca.df$population)

#add fake aus point
pca.df[,1]<-as.numeric(pca.df[,1])
pca.df[,2]<-as.numeric(pca.df[,2])
pcaplot<-ggplot(pca.df,aes(x=X1,y=X2))+science_theme+
  geom_point(aes(colour=population),size=0.7)+
  scale_shape_manual(values=c('circle','square','triangle','plus'))+
  theme(panel.grid=element_blank())+xlim(c(-0.1,0.15))+ylim(c(-0.15,0.1))+
  theme_minimal()+theme(legend.position='none')+
  theme(panel.grid=element_blank(),axis.text=element_blank())+
  xlab("PC1")+ylab("PC2")+labs(tag="A")

#######
# TREE
#######

tree<-read.tree("./fasttree/fasttree_tree_file.txt")
pops<-read.table("./fasttree/pops.txt")

#change pop labels
pops$V2<-factor(pops$V2)
levels(pops$V2)<-c("Cairns","Hawaii.CL","Hawaii.UH","Kauai.AS","Kauai.CG","Kauai.PK",
                   "Kauai.VL","Kauai.WC","Mission","Oahu.AC","Oahu.BYU","Oahu.CC","Oahu.KP","Outgroup")
pops$island<-"commodus"
pops[grep("Kauai",pops$V2),]$island<-"Kauai"
pops[grep("Oahu",pops$V2),]$island<-"Oahu"
pops[grep("Hawaii",pops$V2),]$island<-"Hawaii"
pops[grep("Mission",pops$V2),]$island<-"Australia"
pops[grep("Cairns",pops$V2),]$island<-"Australia"

drop<-tree$tip.label[!tree$tip.label %in% pops$V1]
tree<-drop.tip(tree,drop)
pops[pops$island=="Australia" | pops$island=="commodus_outgroup",]$V2<-"outgroup/Australia"

#root tree by Australian mainland oceanicus sample
tree<-root(tree,outgroup = 'A4')
p <- ggtree(tree,size=0.5,aes(colour=V2),linewidth=0.5) %<+% pops +
  geom_tippoint(aes(color = V2), size = 0.4) +  
  theme_tree2() +  science_theme +
  theme(legend.title = element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        panel.background=element_rect(colour='white'),panel.border=element_blank(),
        legend.position = 'none')+
  labs(tag="B")

p<-flip(p, 377, 313)
print(p)


#######
# FST
#######

fst<-read.csv("./pairwise_fst_pops/fst.csv",row.names = 1)

#create and plot numeric matrix
fst<-apply(fst,1,as.numeric)
rownames(fst)<-colnames(fst)
#removal diagonal
diag(fst)<-NA
fst.heatmap<-pheatmap(fst,cluster_cols = FALSE,cluster_rows = FALSE,gaps_row = c(2,7),
                      gaps_col = c(2,7),show_rownames = TRUE,show_colnames = TRUE,
                      border_color = '#555555',display_numbers = FALSE,fontsize = 6)
#convert to ggplot object for multi-panel plot
fst.heatmap<-as.ggplot(fst.heatmap)+labs(tag="C")

#######
# PHLASH
#######

q.0.25<-function(x){quantile(x,0.025)}
q.0.975<-function(x){quantile(x,0.975)}

phlash<-read.csv("./phlash/Nes_phlash_OahuBYU_chr5_9_14.csv",h=F)
phlash<-data.frame(t(phlash))
oahuBYU.t<-read.csv("./phlash/Nes_phlash_OahuBYU_chr5_9_14_timepoints.csv",h=F)
phlash$pop<-"OahuBYU"
phlash$T<-oahuBYU.t$V1

#each column is a different iteration, so take the median for each sampled point (row)
phlash[,c(1:500)]<-data.frame(apply(phlash[,c(1:500)],2,as.numeric))
phlash$median<-(apply(phlash[,c(1:500)],1,median))
phlash$q5<-apply(phlash[,c(1:500)],1,q.0.25)
phlash$q95<-apply(phlash[,c(1:500)],1,q.0.975)

phlash.plot<-ggplot(phlash,aes(x=(T),y=median))+
  science_theme+theme(panel.grid=element_blank(),axis.ticks=element_line())+
  geom_ribbon(aes(ymin = (q5),ymax=(q95)),alpha=0.15)+
  geom_vline(xintercept=827*3.5,colour='black',linetype='dashed')+
  geom_vline(xintercept=116,colour='darkred',linetype='dotted',linewidth=0.75)+
  geom_line(size=1)+
  scale_y_log10()+scale_x_log10()+
  xlab("Generation")+annotation_logticks()+ylab("Effective population size")+
  labs(tag="D")

#######
# GONE2
#######

setwd("./gone2")

gone2<-rbind(
  read.table("Hilo.CL_r1.5_full_GONE2_Ne",h=T) %>% mutate(pop="Hilo.CL") %>% mutate(island="Hawaii"),
  read.table("Hilo.UH_r1.5_full_GONE2_Ne",h=T) %>% mutate(pop="Hilo.UH") %>% mutate(island="Hawaii"),
  read.table("Kauai.AS_r1.5_full_GONE2_Ne",h=T) %>% mutate(pop="Kauai.AS") %>% mutate(island="Kauai"),
  read.table("Kauai.VL_r1.5_full_GONE2_Ne",h=T) %>% mutate(pop="Kauai.VL") %>% mutate(island="Kauai"),
  read.table("Kauai.CG_r1.5_full_GONE2_Ne",h=T) %>% mutate(pop="Kauai.CG") %>% mutate(island="Kauai"),
  read.table("Kauai.PK_r1.5_full_GONE2_Ne",h=T) %>% mutate(pop="Kauai.PK") %>% mutate(island="Kauai"),
  read.table("Oahu.AC_r1.5_full_GONE2_Ne",h=T) %>% mutate(pop="Oahu.AC") %>% mutate(island="Oahu"),
  read.table("Oahu.BYU_r1.5_full_GONE2_Ne",h=T) %>% mutate(pop="Oahu.BYU") %>% mutate(island="Oahu"),
  read.table("Oahu.CC_r1.5_full_GONE2_Ne",h=T) %>% mutate(pop="Oahu.CC") %>% mutate(island="Oahu"),
  read.table("Oahu.KP_r1.5_full_GONE2_Ne",h=T) %>% mutate(pop="Oahu.KP") %>% mutate(island="Oahu")
)


g.ne.full<-ggplot(gone2,aes(x=Generation,y=(Ne_diploids),colour=pop))+science_theme+
  geom_rect(xmin=42,xmax=66.5,ymin=min(log10(gone2$Ne_diploids)),
            ymax=max(log10(gone2$Ne_diploids)),fill='#eeeeee',colour='#eeeeee')+
  geom_vline(xintercept=115.5,colour='darkred',linetype='dotted',linewidth=0.75)+
  geom_line(linewidth=0.66,aes(fill=pop),alpha=1,linetype='solid')+
  stat_summary(geom='line',fun='median',colour='black',linewidth=1.25,show.legend = FALSE)+
  ylab("Ne")+guides(color = guide_legend(override.aes = list(linewidth = 2)))+
  #facet_wrap(.~pop)+
  theme(panel.grid=element_blank(),panel.background=element_rect(fill='white',colour='black'),
        legend.position='right',strip.background=element_rect(fill="#eeeeee"))+
  xlab("generation")+scale_y_log10()+annotation_logticks(sides = 'l')+ylab("Effective population size")+
  labs(tag="E")



library(patchwork)

ggsave("fig3_p1.svg",height=3,width=7.25,
       plot=(pcaplot+p+fst.heatmap+plot_layout(widths=c(0.85,0.66,1))))

ggsave("fig3_p2.svg",height=3,width=7.25,
       plot=(phlash.plot+g.ne.full))


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
#   [1] pheatmap_1.0.12  plyr_1.8.9       ggbeeswarm_0.7.2 lubridate_1.9.4  forcats_1.0.0    stringr_1.5.1   
# [7] purrr_1.0.4      readr_2.1.5      tibble_3.2.1     tidyverse_2.0.0  tidyr_1.3.1      ggplotify_0.1.2 
# [13] pcadapt_4.4.0    patchwork_1.3.0  ggExtra_0.10.1   ape_5.8-1        dplyr_1.1.4      ggplot2_4.0.2   
# [19] ggtree_3.12.0   
# 
# loaded via a namespace (and not attached):
#   [1] yulab.utils_0.2.0  generics_0.1.3     stringi_1.8.4      lattice_0.22-6     hms_1.1.3         
# [6] digest_0.6.37      magrittr_2.0.3     timechange_0.3.0   grid_4.4.2         RColorBrewer_1.1-3
# [11] fastmap_1.2.0      jsonlite_1.8.9     promises_1.3.2     aplot_0.2.4        scales_1.4.0      
# [16] lazyeval_0.2.2     cli_3.6.3          shiny_1.10.0       rlang_1.1.5        tidytree_0.4.6    
# [21] withr_3.0.2        tools_4.4.2        parallel_4.4.2     tzdb_0.4.0         httpuv_1.6.15     
# [26] vctrs_0.6.5        R6_2.5.1           mime_0.12          gridGraphics_0.5-1 lifecycle_1.0.4   
# [31] fs_1.6.5           ggfun_0.1.8        vipor_0.4.7        miniUI_0.1.1.1     treeio_1.28.0     
# [36] beeswarm_0.4.0     pkgconfig_2.0.3    pillar_1.10.1      later_1.4.1        gtable_0.3.6      
# [41] glue_1.8.0         Rcpp_1.1.0         tidyselect_1.2.1   rstudioapi_0.17.1  farver_2.1.2      
# [46] xtable_1.8-4       htmltools_0.5.8.1  nlme_3.1-167       compiler_4.4.2     S7_0.2.1  