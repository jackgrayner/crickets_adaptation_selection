library(ggplot2)
library(tidyverse)
library(reshape2)

pca.1<-data.frame()
for (chr in c(1:14)){
  pca <- read.table(paste0("./winpca_results/merged_scaffold_",chr,".pc_1.tsv"),h=T) %>% melt(id.var=c("pos"))
  pca$chr<-chr
  pca.1<-rbind(pca.1,pca)
}

# pca.2<-data.frame()
# for (chr in c(2:14)){
#   pca <- read.table(paste0("merged_scaffold_",chr,".pc_2.tsv"),h=T) %>% melt(id.var=c("pos"))
#                     pca$chr<-chr
#                     pca.2<-rbind(pca.2,pca)
# }


pca.1$pop="Hawaii"
pca.1[grep("ERR",pca.1$variable),]$pop="Australia"
pca.1$col<-paste(pca.1$pop,pca.1$chr %% 2 ==0)
summary(factor(pca.1$col))

pca.1$alpha=0.5
pca.1[pca.1$chr %% 2 ==0,]$alpha=1

g.pc1.space<-ggplot(pca.1,aes(x=pos,y=value,colour=col))+theme_bw()+
  theme(panel.grid=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),legend.position='none')+
  geom_line(aes(group=variable),linewidth=0.025,alpha=0.75)+
  facet_grid(pop~chr,scale='free',space='free')+
  theme(panel.spacing = unit(0.1, "lines"),
        strip.text = element_text(margin = margin(0.05,0.05,0.05,0.05, "cm")),
        #panel.border=element_rect(fill='#fff',colour="#888"),
        strip.background=element_rect(fill='#eeeeee',colour="#333"))+
  #scale_colour_viridis(discrete = TRUE,end = 0.7)+
  scale_colour_manual(values=c("#8D918B","#8D918B","#99c9a0","#99c9a0"))+
  scale_x_continuous(expand = c(0, 0))+
  ylab("PC1")

ggsave("pca1.png",plot=g.pc1.space+labs(tag=" "),dpi=600,height=2.25,width=7.86)


#chrs of interest
min(pca.1[pca.1$chr==10,]$pos)
g.chrs.oi<-ggplot(pca.1[pca.1$chr %in% c(1,4,6,7,10),],aes(x=pos/1e+06,y=value,colour=col))+theme_bw()+
  theme(panel.grid=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),legend.position='none')+
  geom_line(aes(group=variable),linewidth=0.05,alpha=1)+
  facet_grid(pop~chr,scale='free',space='free')+
  theme(axis.text.x=element_text(size=6),
        strip.text = element_text(margin = margin(0.05,0.05,0.05,0.05, "cm")),
        #panel.border=element_rect(fill='#fff',colour="#888"),
        strip.background=element_rect(fill='#eeeeee',colour="#333"))+
  #scale_colour_viridis(discrete = TRUE,end = 0.7)+
  scale_colour_manual(values=c("#8D918B","#8D918B","#99c9a0","#99c9a0"))+
  scale_x_continuous(breaks=seq(0,350,50))+
  ylab("PC1")

ggsave("pca1_chrs_of_interest.svg",plot=g.chrs.oi,dpi=600,height=2.15,width=4.5)


#add selected regions
#run fig 4 code first
flysel.labels$chr.x<-flysel.labels$chr
g.pc1<-ggplot(pca1.1[pca1.1$pop=="HI",],aes(x=pos/1e+06,y=value,colour=scaled.mean))+theme_bw()+
  geom_vline(data=flysel.labels,aes(xintercept=ps/1e+06),colour='#b31e2c',linewidth=0.65)+
  theme(panel.grid=element_blank(),axis.title.x=element_blank(),legend.position='none',strip.text.y = element_text(face = 'bold'))+
  geom_line(aes(group=variable.x),linewidth=0.1,alpha=0.8)+
  facet_wrap(.~chr.x,scale='free')+
  theme(panel.spacing = unit(0, "lines"))+
  scale_colour_viridis(option = "D")+
  scale_x_continuous(expand = c(0, 0))+
  ylab("PC1")

ggsave("pca1_winpc_sel_snps.png",plot=g.pc1,dpi=600,height=8,width=9)


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
#   [1] reshape2_1.4.4  lubridate_1.9.4 forcats_1.0.0   stringr_1.5.1   dplyr_1.1.4     purrr_1.0.4    
# [7] readr_2.1.5     tidyr_1.3.1     tibble_3.2.1    tidyverse_2.0.0 ggplot2_4.0.2  
# 
# loaded via a namespace (and not attached):
#   [1] gtable_0.3.6       compiler_4.4.2     tidyselect_1.2.1   Rcpp_1.1.0         scales_1.4.0      
# [6] R6_2.5.1           plyr_1.8.9         generics_0.1.3     pillar_1.10.1      RColorBrewer_1.1-3
# [11] tzdb_0.4.0         rlang_1.1.5        stringi_1.8.4      S7_0.2.1           timechange_0.3.0  
# [16] cli_3.6.3          withr_3.0.2        magrittr_2.0.3     grid_4.4.2         rstudioapi_0.17.1 
# [21] hms_1.1.3          lifecycle_1.0.4    vctrs_0.6.5        glue_1.8.0         farver_2.1.2      
# [26] tools_4.4.2        pkgconfig_2.0.3  