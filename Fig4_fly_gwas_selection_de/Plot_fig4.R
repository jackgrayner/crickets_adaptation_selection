library(dplyr)
library(ggplot2)
library(patchwork)
library(ggrepel)
library(topGO)
library(GenomicRanges)
library(ggbeeswarm)
library(ggplot2)
library(viridis)
library(lme4)
library(car)
library(rptR)
library(dplyr)


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
# fly plots
######

flies<-read.csv('./Fly_selection_rates/flytrapping.csv',h=T)
flies<-flies[flies$rain %in% c("dry","light rain"),]#get rid of trials during heavy rain

sum(flies[flies$playback=="song",]$flies_caught)
sum(flies[flies$playback=="control",]$flies_caught)

#relabel sites
flies$site<-gsub("K\\.","Kauai.",flies$site)
flies$site<-gsub("O\\.","Oahu.",flies$site)
flies$site<-gsub("H\\.","Hawaii.",flies$site)

#classify sites according to proportion of singing crickets
flies$wtnwpres<-"no_song"
flies[flies$site %in% c("Kauai.CG","Kauai.HC","Oahu.CC","Oahu.KP"),]$wtnwpres<-"low_song"
flies[flies$site %in% c("Hawaii.CL","Kauai.PK","Hawaii.UH","Oahu.BYU"),]$wtnwpres<-"high_song"
flies$wtnwpres<-factor(flies$wtnwpres,levels=c("no_song","low_song","high_song"))

#create binary variable (0 = 0 flies, 1 = >0 flies)
flies$fly.yn<-0
flies[flies$flies_caught>0,]$fly.yn<-1

flies$island<-factor(flies$island,levels=c("Kauai","Oahu","Hawaii"))
g.time<-ggplot(flies,aes(x=time_postsunset,y=flies_caught))+science_theme+
  geom_rect(fill='#fcf3f2',xmin=0,xmax=120,ymin=0,ymax=1.25)+coord_cartesian(ylim=c(0,1.25))+
  theme(legend.position='none')+
  geom_smooth(method='loess',span=1,aes(group=date,colour=island),se=FALSE,linewidth=0.35)+
  geom_smooth(method='loess',span=1,colour='black',alpha=1,fill='#cccccc',size=1.5)+
  ylab("Flies caught per trap")+xlab('Mins. post-sunset')+
  scale_colour_manual(values=c("#00c08b","#c77cff","#f8766d"))

#now estimate relative fitness of silent males

#1. take mean flies attracted per song trap
flies.sum.night<-flies %>% filter(playback=='song' & time_postsunset<120) %>% 
  group_by(site,island,date) %>% 
  dplyr::summarise(meanflies=mean(flies_caught),meanfly.yn=mean(fly.yn))
unique(flies.sum.night$site)

#calculate repeatability of mean likelihood of attracting a fly (Y/N) across sites 
#rep1<-rpt(formula=(log(meanfly.yn+1)~(1|site)),grname="site",data=data.frame(flies.sum.night),datatype = "Gaussian")
#qqnorm(resid(rep1$mod))#resids look reasonably normally distributed
#qqline(resid(rep1$mod))
#shapiro.test(resid(rep1$mod))#as above
#plot(rep1$mod)#could be heteroscedasticity
#leveneTest(log(meanfly.yn+1) ~ factor(site), data = data.frame(flies.sum.night))#P = 0.10. Probably borderline 
#plot(rep1)

#test for island/site diffs in likelihood of attracting fly
glm2<-lm(log(meanfly.yn+1)~island/site+season,data=data.frame(flies.sum.night))
anova(glm2)
car::Anova(glm2)
plot(glm2)

#that's for 10 mins (Pten). crickets sing for 7.63 mins (Pfly)
#'expected arrivals' per 7.63 minutes is
Nten = flies.sum.night$meanflies
Nfly = Nten*0.76

#probability of  being infected per night is...
#(only 61% of attacks end in infestation)
Pinf = 1 - (1 - 0.61)^(Nfly)

est_lifespan = 1/Pinf + 6 #add 6 day lag
est_lifespan[est_lifespan>29]<-29 #add 29 day expected lifespan baseline
flies.sum.night$lifeexpectency<-est_lifespan
mean(est_lifespan)
sd(est_lifespan)

NS_relfitness_fw<-29/(flies.sum.night$lifeexpectency)
mean(NS_relfitness_fw)
sd(NS_relfitness_fw)
SS_relfitness_fw<-(1-0.391)#Sexual selection estimate taken from Tanner 2019
flies.sum.night$fw.rel.fitness<-NS_relfitness_fw*SS_relfitness_fw # = net rel. fitness of Fw males
mean(flies.sum.night$fw.rel.fitness)
sd(flies.sum.night$fw.rel.fitness)
min((flies.sum.night$fw.rel.fitness))
max((flies.sum.night$fw.rel.fitness))

flies.sum.night$island<-factor(flies.sum.night$island,levels=c("Kauai","Oahu","Hawaii"))
g.silent<-ggplot(flies.sum.night,aes(x=site,y=log2(fw.rel.fitness),colour=island,group=date))+
  science_theme+
  geom_hline(yintercept=0,linetype='dashed')+
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1),
                   axis.title.x=element_blank(),legend.position = 'none',
                   strip.background=element_rect(fill='#eeeeee'))+
  facet_grid(.~island,scales='free',space='free')+geom_quasirandom(alpha=0.75,width=0.1)+
  scale_y_continuous(breaks=seq(-1,2,0.2))+
  ylab("log2 rel. fitness of silent males")+#scale_y_continuous(transform = 'log2')+annotation_logticks(sides = 'l')+
  stat_summary(aes(group=site,fill=island),colour='#555555',shape=21,linewidth=0.5,
               fun.data=mean_sdl,size=0.5,
               fun.args = list(mult = 1))+
  #geom_boxplot(aes(group=site),alpha=0)+
  scale_colour_manual(values=c("#00c08b","#c77cff","#f8766d"))+
  scale_fill_manual(values=c("#00c08b","#c77cff","#f8766d"))#+
#scale_y_continuous(breaks=c(-0.8,1.2,0.2))


#######
#read/filter genome annotation, add functional annotations
######

#read genome annotation, reformat some columns, subset genes
tocgff<-read.table("./GWAS_selection_analysis/TOC.asm.scaffold.gene.gff3") %>% 
  mutate(V1=as.integer(gsub("scaffold_","",V1)),
   #      gene=gsub("ID=.*;.*Name=","",V9)) %>%
  #filter(V3=="gene")
  gene=gsub("ID=exon.*;.*Parent=","",V9)) %>%
  filter(V3=="exon")
tocgff$gene<-gsub("\\.t1*","",tocgff$gene)


#read functional annotations
gene_annos<-read.table("./GWAS_selection_analysis/TOC.asm.scaffold.gene.SWISSPROT.blastp.top_hit.txt")
colnames(gene_annos)<-c("gene","score","anno")
gene_annos$gene<-gsub("\\.t.*","",gene_annos$gene)#remove transcript identifier

#keep top annotation per gene
gene_annos<-gene_annos[order(gene_annos$score,decreasing=TRUE),] %>% filter(!duplicated(gene_annos))

#remove accession info and spp ID from annotations for readability, then merge genome and functional annotation data frames
gene_annos$anno<-gsub(".*\\|","",gene_annos$anno)
#gene_annos$anno<-gsub("_.*","",gene_annos$anno)
tocgff.anno<-merge(tocgff,gene_annos,by='gene')
tocgff.anno[tocgff.anno$anno=="DSX_DROME",]
tocgff.anno[tocgff.anno$anno=="WLS_DROPS",]

#create granges object used to identify overlapping regions (i.e. mapping variants to nearby genes)
genes<-GRanges(seqnames=tocgff.anno$V1,#chr/scaffolds
  ranges=IRanges(start=tocgff.anno$V4, end=tocgff.anno$V5))#start and end of gene feature
mcols(genes)$gene <- tocgff.anno$gene#add gene ID (e.g., TOC.11T.00352)
mcols(genes)$anno <- tocgff.anno$anno#add annotation (eg., UBCD1)

window=50000 #50kb window used to identify nearby genes
#find nearby genes for sig. variants, output sig. vars with genes & annotations
annotation_function<-function(df){ 
  #make empty columns
  df$gene<-NA#to store gene IDs for nearby genes
  df$anno<-NA#to store annotations for nearby genes
  df$dist<-NA#to store distance to the nearby gene (0 if overlapping)
  df$within<-FALSE#whether the variant is within a gene feature
  
  #convert input into granges object
  variants<-GRanges(seqnames=df$chr,ranges=IRanges(start=df$ps,end=df$ps))
  #find overlaps between query/subject (i.e. variants within genes)
  variants_within<-findOverlaps(variants, genes)
  
  #store within (TRUE) and distance (0) for overlapping ranges
  df$within[queryHits(variants_within)]<-TRUE
  df$dist[queryHits(variants_within)]<-0
  
  #find nearby genes
  dist.nearest<-distanceToNearest(variants,genes,ignore.strand=TRUE)
  nearest_genes<-genes[subjectHits(dist.nearest)]
  distances<-mcols(dist.nearest)$distance
  #save results to variants df
  df$gene<-mcols(nearest_genes)$gene
  df$anno<-mcols(nearest_genes)$anno
  df$dist<-distances
  #note whether dist is < 50Kb
  df$within50kb<-df$dist<window
  return(df)
}


#######
#curly-wing GWAS - plot (only Chr2)
######

#read subset of association test results, keep only Chr2
All.Cw<-read.table("./GWAS_selection_analysis/All_Cw_lmm.assoc_P0.1.txt",h=T) %>% filter(chr==2)
topcw<-All.Cw[order(All.Cw$p_lrt),][1,]$ps

#annotate sig. SNPs
All.Cw$sig<-All.Cw$padj<0.05
All.Cw<-annotation_function(All.Cw)

cw.top<-head(All.Cw[order(All.Cw$p_lrt),],n=10)
cw.top$chr<-paste0("scaffold_",cw.top$chr)

#create labelled gene set.
cw.labels<-All.Cw[All.Cw$sig  & All.Cw$within50kb,]
cw.labels<-cw.labels[order(cw.labels$p_lrt),]
cw.labels<-cw.labels[!duplicated(cw.labels$anno),]

#write.csv(file="~/Documents/StA/popgen/gwas/Rayner2025_Cw_top100.csv",head(All.Cw[order(All.Cw$p_lrt),],n=100),quote=FALSE,row.names = FALSE)

#plot Cw association across Chr2
All.Cw.plot<-ggplot(All.Cw,aes(x=ps/1e+6,y=(-log10(p_lrt)),colour=sig))+
  xlab("Chr2 pos. Mb")+science_theme+theme(axis.ticks=element_line())+
  geom_point(size=0.75,alpha=1)+
  #plot gene annotations for SNPs that are significant and within 50kb of genes
  geom_text_repel(data=cw.labels[cw.labels$padj<0.00001,],
                  aes(x=ps/1e+6,y=(-log10(p_lrt)),label=gsub("_.*","",anno)),
                  min.segment.length=0.0001,segment.colour='#aaaaaa',size=2,colour='black',
                  nudge_x=0,nudge_y=1.5,fontface="italic",alpha=1,segment.size=0.25)+
  scale_colour_manual(values=c("#d4cdcd","#8f3131"))+
  ggtitle("Curly-wing")+ylab("-log10(P)")+
  scale_y_continuous(expand=c(0,0))

rm(list="All.Cw")#remove from environment to save memory

#######
#flatwing GWAS - plot (only Chr1)
######

#read subset of association test results for each pop, keep only Chr1
kauai.fw<-read.table("./GWAS_selection_analysis/Kauai_Fw.assoc_P0.1.txt",h=T) %>% filter(chr=="scaffold_1") %>% mutate(chr=1,sig=Padj<0.05)
oahu.fw<-read.table("./GWAS_selection_analysis/Oahu_Fw.assoc_P0.1.txt",h=T) %>% filter(chr=="scaffold_1") %>% mutate(chr=1,sig=Padj<0.05)
#hilo.fw<-read.table("Hawaii_Fw.assoc_P0.1.txt",h=T) %>% filter(chr=="scaffold_1") #excl. as too few Fw samples

#store positions of top variants for plotting
topfw.kauai<-kauai.fw[order(kauai.fw$p_lrt),][1,]$ps
topfw.oahu<-oahu.fw[order(oahu.fw$p_lrt),][1,]$ps

kfw.top<-head(kauai.fw[order(kauai.fw$p_lrt),],n=10)
kfw.top$chr<-paste0("scaffold_",kfw.top$chr)
write.table(kfw.top[,c("chr","ps")],quote=FALSE,row.names=FALSE,col.names=FALSE,
            file="kfw.top.tsv",sep="\t")
ofw.top<-head(kauai.fw[order(oahu.fw$p_lrt),],n=10)
ofw.top$chr<-paste0("scaffold_",ofw.top$chr)
write.table(ofw.top[,c("chr","ps")],quote=FALSE,row.names=FALSE,col.names=FALSE,
            file="ofw.top.tsv",sep="\t")

#add annotations
kauai.fw<-annotation_function(kauai.fw)
oahu.fw<-annotation_function(oahu.fw)

#create labelled gene set. remove duplicated genes independently 
kauai.fw.labels<-kauai.fw[kauai.fw$sig  & kauai.fw$within50kb,]
kauai.fw.labels<-kauai.fw.labels[order(kauai.fw.labels$p_lrt),]
kauai.fw.labels<-kauai.fw.labels[!duplicated(kauai.fw.labels$anno),]
oahu.fw.labels<-oahu.fw[oahu.fw$sig  & oahu.fw$within50kb,]
oahu.fw.labels<-oahu.fw.labels[order(oahu.fw.labels$p_lrt),]
oahu.fw.labels<-oahu.fw.labels[!duplicated(gsub("_.*","",oahu.fw.labels$anno)),]
fw.labels<-rbind(kauai.fw.labels,oahu.fw.labels)

#concatenate
All.Fw<-rbind(kauai.fw,oahu.fw)
All.Fw[All.Fw$rs %in% kauai.fw[kauai.fw$sig,]$rs & All.Fw$rs %in% oahu.fw[oahu.fw$sig,]$rs,]$sig<-"Both"
mean(All.Fw[All.Fw$rs %in% kauai.fw[kauai.fw$sig,]$rs & All.Fw$rs %in% oahu.fw[oahu.fw$sig,]$rs,]$ps)
rm(list=c("kauai.fw","oahu.fw"))#save memory

#define whether SNP is significant in either island
All.Fw$sig<-NA
All.Fw[All.Fw$Padj<0.05 & All.Fw$Island=="Kauai",]$sig<-"Kauai"
All.Fw[All.Fw$Padj<0.05 & All.Fw$Island=="Oahu",]$sig<-"Oahu"

nrow(All.Fw[All.Fw$sig=="Both",])
fw.both<-All.Fw[All.Fw$sig=="Both",]
fw.both<-fw.both[!duplicated(fw.both$rs),]
fw.both$chr<-paste0("scaffold_",fw.both$chr)

#plot Fw associations across chr1, split by island
All.Fw.plot<-ggplot(All.Fw,aes(x=ps/1e+6,y=(-log10(p_lrt)),colour=sig))+facet_grid(Island~.,scales='free')+
  xlab("Chr1 pos. Mb")+science_theme+theme(strip.background=element_rect(fill='#eeeeee'))+
  geom_vline(xintercept=253.6,linetype='dashed',linewidth=0.25,colour='black')+#add dsx location
  geom_point(size=0.75,alpha=1)+
  #plot gene annotations for SNPs that are significant and within 50kb of genes
  geom_text_repel(data=fw.labels,segment.size=0.25,
                  aes(x=ps/1e+6,y=(-log10(p_lrt)),label=gsub("_.*","",anno)),colour='black',max.overlaps=10,
                  min.segment.length=0.0001,segment.colour='#aaaaaa',size=2,nudge_x=0,nudge_y=1.5,fontface="italic",alpha=1)+
  scale_colour_manual(values=c("#00c08b","#c77cff","red"),na.value="#d4cdcd")+
  ggtitle("Flatwing")+ylab("-log10(P)")+
  scale_y_continuous(expand=c(0,0))#+coord_cartesian(xlim=(c(253,255)),ylim=c(0,25))

rm(list="All.Fw")

#######
#selection results (LFMM of genotype ~ fly attack rate)
######
flysel<-read.csv("./GWAS_selection_analysis/LFMM_Pinf_results.csv",h=T) %>% mutate(p_lrt=P,chr=as.integer(gsub("scaffold_","",chr))) %>% filter(!is.na(Padj))
flysel$sig<-flysel$Padj<0.05
flysel$ps1<-round(flysel$ps/2.5e+06)*2.5
flysel$chr1.ps1<-paste(flysel$chr,flysel$ps1,sep=".")
length(unique(flysel$chr1.ps1))

#add gene info 
flysel<-annotation_function(flysel)

flysel.labels<-flysel[flysel$sig,]
flysel.labels<-flysel.labels[order(flysel.labels$p_lrt),]
flysel.labels<-flysel.labels[!duplicated(flysel.labels$anno),]
flysel.labels$ps1<-round(flysel.labels$ps/10e+06)*50
flysel.labels$chr1.ps1<-paste(flysel.labels$chr,flysel.labels$ps1,sep=".")
length(unique(flysel.labels$chr1.ps1))
flysel.labels<-flysel.labels[!duplicated(flysel.labels$chr1.ps1),]

#plot LFMM results across genome
fw<-median(flysel[flysel$gene =="TOC.1T.01919",]$ps)
cw<-median(flysel[flysel$gene %in% c("TOC.2T.01818","TOC.2T.01817"),]$ps)

fly.sel.plot<-ggplot(flysel,aes(x=ps,y=(-log10(P))))+
  geom_vline(data=data.frame(chr=1),aes(xintercept=fw),linetype='dashed',linewidth=0.35,colour='#aaaaaa')+
  geom_vline(data=data.frame(chr=2),aes(xintercept=cw),linetype='dashed',linewidth=0.35,colour='#aaaaaa')+
  #plot gene annotations for SNPs that are outliers and within 50kb of genes
  science_theme+
  facet_grid(.~chr,space='free',scales='free',switch='both')+
  theme(axis.text.x=element_blank(),axis.title.x=element_blank(),
        panel.spacing.x=unit(0.05, "lines"),panel.border=element_blank(),
        panel.background=element_rect(colour='white',fill='white'),
        strip.background=element_rect(colour='white',fill='white'),
        axis.ticks.x=element_blank())+
  geom_point(aes(colour=sig),size=0.75,alpha=1)+
  geom_text_repel(data=flysel.labels[flysel.labels$sig,],aes(x=ps,y=(-log10(P)),label=gsub("_.*","",anno)),
                  min.segment.length=0.0001,segment.colour='#aaaaaa',size=2,nudge_x=0,nudge_y=1.5,
                  fontface="italic",alpha=1,max.overlaps=1,segment.size=0.25)+
  scale_colour_manual(values=c("#d4cdcd","#f26d79",'#b31e2c'))+
  ggtitle("Fly selection")+ylab("-log10(P)")+scale_y_continuous(expand=c(0,0))
  

######
# DE plots
######

library(DESeq2)
library(topGO)
library(edgeR)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(stringr)
library(variancePartition)

pheno<-read.csv('./Diff_expression/Sikkink_pheno.csv')
rownames(pheno)<-pheno$ID

cts<-read.csv('./Diff_expression/Sikkink_gene_count_matrix.csv',h=T,row.names="gene_id")
summary(colnames(cts) %in% rownames(pheno))
cts<-cts[,order(colnames(cts))]
pheno<-pheno[order(rownames(pheno)),]
summary(colnames(cts) == rownames(pheno))

#DE analysis
dds <- DESeq(DESeqDataSetFromMatrix(countData = cts,
                                    colData = pheno,
                                    design= ~ Pop+Trt))
resultsNames(dds)
summary(results(dds, name="Trt_D4_vs_C",alpha=0.05))
summary(results(dds, name="Trt_D7_vs_C",alpha=0.01))
res.d4<-data.frame(results(dds, name="Trt_D4_vs_C",alpha=0.05))
res.d4$gene<-rownames(res.d4)
res.d7<-data.frame(results(dds, name="Trt_D7_vs_C",alpha=0.05))
res.d7$gene<-rownames(res.d7)

gene_list <- res.d4$padj
names(gene_list) <- rownames(res.d4)
gene_list <- gene_list[!is.na(gene_list)]
selection_function <- function(x) {
  return(x < 0.05)
}
gene2GO <- list()
for(i in 1:nrow(gaf)) {
  gene_symbol <- gaf[i, 1]  # First column is gene symbol
  go_terms <- gaf[i, 2]     # Second column is GO terms
  
  # Skip if GO terms are empty or NA
  if(!is.na(go_terms) && go_terms != "") {
    # Split GO terms by semicolon and clean whitespace
    go_list <- trimws(unlist(strsplit(go_terms, ";")))
    # Remove empty strings
    go_list <- go_list[go_list != ""]
    gene2GO[[gene_symbol]] <- go_list
  }
}

GOdata_BP <-  new("topGOdata",
                  ontology = "BP",
                  allGenes = gene_list,
                  geneSel = selection_function,
                  annot = annFUN.gene2GO,
                  gene2GO = gene2GO)
GOdata_MF <- new("topGOdata",
                 ontology = "MF",
                 allGenes = gene_list,
                 geneSel = selection_function,
                 annot = annFUN.gene2GO,
                 gene2GO = gene2GO)

test_BP <- runTest(GOdata_BP, algorithm = "weight01", statistic = "fisher")
test_MF <- runTest(GOdata_MF, algorithm = "weight01", statistic = "fisher")
results_BP <- GenTable(GOdata_BP, Fisher = test_BP, topNodes = 100)
results_MF <- GenTable(GOdata_MF, Fisher = test_MF, topNodes = 100)
results_BP$Category <- "Biological Process"
results_MF$Category <- "Molecular Function"
all_results <- rbind(results_BP, results_MF)
all_results$Fisher<-as.numeric(all_results$Fisher)
all_results<-all_results[all_results$Significant>1,]
all_results$Padj<-p.adjust(all_results$Fisher)
all_results$fold<-all_results$Significant/all_results$Annotated

#plot
g.volcano<-ggplot(res.d4[!is.na(res.d4$padj),],aes(x=log2FoldChange,y=-log10(pvalue)))+
  science_theme+theme(legend.position='none')+
  geom_point(size=0.85,alpha=0.75,aes(colour=padj<0.05))+
  geom_vline(xintercept=0,linetype='dashed')+
  xlim(c(-12,12))+ylim(c(0,42))+
  scale_colour_manual(values=c("#aaaaaa","#C74955"))+
  scale_fill_manual(values=c("#659dd4","#d59ef7"))

res.d4.go<-all_results %>%
  mutate(set = case_when(
    GO.ID %in% c("GO:0004252","GO:0005178") ~ "serps.integ",
    GO.ID %in% "GO:0035001" ~ "tracheas",
  ))

g.go<-ggplot(all_results[all_results$Fisher<3.5e-5,],aes(y=Term,x=-log10(Fisher)))+
  science_theme+theme(legend.position='none',axis.title.y=element_blank(),axis.line.x=element_line())+
  geom_col(colour='black')+xlab("Fold enrichment")#+scale_fill_manual(values=c("#659dd4","#d59ef7"))

######
# final figure
######

ggsave('fig4.png',dpi=1200,height=7.8,width=7.25,
       plot=(g.time+labs(tag="A",title="Fly attack rate")+g.silent+labs(tag="B"))/
         (All.Fw.plot+labs(tag="C")+All.Cw.plot+labs(tag="D"))/
         (fly.sel.plot+labs(tag="E"))/
         (g.volcano+labs(tag="F",title="Response to infestation")+g.go+labs(tag="G"))+
         plot_layout(nrow=4,heights=c(5,4,4,5)))


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
#   [1] variancePartition_1.34.0    BiocParallel_1.38.0         stringr_1.5.1              
# [4] RColorBrewer_1.1-3          pheatmap_1.0.12             edgeR_4.2.2                
# [7] limma_3.60.6                DESeq2_1.44.0               SummarizedExperiment_1.34.0
# [10] MatrixGenerics_1.16.0       matrixStats_1.5.0           rptR_0.9.22                
# [13] car_3.1-3                   carData_3.0-5               lme4_1.1-36                
# [16] Matrix_1.7-2                viridis_0.6.5               viridisLite_0.4.2          
# [19] ggbeeswarm_0.7.2            GenomicRanges_1.56.2        GenomeInfoDb_1.40.1        
# [22] topGO_2.56.0                SparseM_1.84-2              GO.db_3.19.1               
# [25] AnnotationDbi_1.66.0        IRanges_2.38.1              S4Vectors_0.42.1           
# [28] Biobase_2.64.0              graph_1.82.0                BiocGenerics_0.50.0        
# [31] ggrepel_0.9.6               patchwork_1.3.0             ggplot2_4.0.2              
# [34] dplyr_1.1.4                
# 
# loaded via a namespace (and not attached):
#   [1] Rdpack_2.6.2            DBI_1.2.3               bitops_1.0-9            gridExtra_2.3          
# [5] rlang_1.1.5             magrittr_2.0.3          compiler_4.4.2          RSQLite_2.3.9          
# [9] reshape2_1.4.4          png_0.1-8               vctrs_0.6.5             pkgconfig_2.0.3        
# [13] crayon_1.5.3            fastmap_1.2.0           backports_1.5.0         XVector_0.44.0         
# [17] caTools_1.18.3          UCSC.utils_1.0.0        nloptr_2.1.1            purrr_1.0.4            
# [21] bit_4.5.0.1             zlibbioc_1.50.0         cachem_1.1.0            jsonlite_1.8.9         
# [25] EnvStats_3.0.0          blob_1.2.4              remaCor_0.0.18          DelayedArray_0.30.1    
# [29] broom_1.0.7             parallel_4.4.2          R6_2.5.1                stringi_1.8.4          
# [33] boot_1.3-31             numDeriv_2016.8-1.1     iterators_1.0.14        Rcpp_1.1.0             
# [37] splines_4.4.2           tidyselect_1.2.1        rstudioapi_0.17.1       abind_1.4-8            
# [41] gplots_3.2.0            codetools_0.2-20        plyr_1.8.9              lmerTest_3.1-3         
# [45] lattice_0.22-6          tibble_3.2.1            withr_3.0.2             KEGGREST_1.44.1        
# [49] S7_0.2.1                Biostrings_2.72.1       pillar_1.10.1           KernSmooth_2.23-26     
# [53] reformulas_0.4.0        generics_0.1.3          RCurl_1.98-1.16         scales_1.4.0           
# [57] aod_1.3.3               minqa_1.2.8             gtools_3.9.5            RhpcBLASctl_0.23-42    
# [61] glue_1.8.0              tools_4.4.2             fANCOVA_0.6-1           locfit_1.5-9.11        
# [65] mvtnorm_1.3-3           grid_4.4.2              tidyr_1.3.1             rbibutils_2.3          
# [69] colorspace_2.1-1        nlme_3.1-167            GenomeInfoDbData_1.2.12 beeswarm_0.4.0         
# [73] vipor_0.4.7             Formula_1.2-5           cli_3.6.3               S4Arrays_1.4.1         
# [77] corpcor_1.6.10          gtable_0.3.6            pbkrtest_0.5.3          SparseArray_1.4.8      
# [81] farver_2.1.2            memoise_2.0.1           lifecycle_1.0.4         httr_1.4.7             
# [85] statmod_1.5.0           bit64_4.6.0-1           MASS_7.3-64 