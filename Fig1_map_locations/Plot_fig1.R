library("ggplot2", "sf")
library("rnaturalearth")
library("rnaturalearthdata")
library("ggspatial")
library(viridis)
library(wesanderson)
library(ape)
library(ggrepel)

science_theme <- theme_minimal(base_size = 6) +
  theme(
    panel.grid=element_blank(),
    text = element_text(family = "Arial"),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7.5),
    strip.text = element_text(size = 8),
    plot.tag = element_text(size = 9, face = "bold")
  )

world <- ne_countries(scale = "medium", returnclass = "sf")

#import data 
locations<-read.csv("cw_mapdata_silentvars.csv")
locations$Silent_Phenotypes<-locations$Phenotypes1
locations$fontface<-'bold'
locations[locations$singing_remain=="N",]$fontface<-'italic'
Map_Populations_crickets<-ggplot(data=world)+science_theme+
  geom_sf(fill='antiquewhite',colour="black")+
  geom_point(data=locations,aes(x=longitude,y=latitude,fill=name1),shape=21,size=1.6,alpha=0.8,colour='black')+
  geom_text_repel(data=locations,
   aes(x=longitude,y=latitude,label=name1,colour=name1),
    size=2.8,force=50,show.legend=FALSE,force_pull=0.1,segment.size = 0.25,min.segment.length = 0.1)+
  annotation_scale(location='bl',width_hint=0.5)+
  annotation_north_arrow(location='bl',which_north='true',pad_x=unit(0.75,"in"),
  pad_y=unit(0.5,"in"),style=north_arrow_fancy_orienteering)+
  coord_sf(xlim=c(-160.580614,-154.278239),ylim=c(18.268486,22.798264),expand=FALSE)+
  xlab('Longitude')+ylab('Latitude')+
  theme(panel.grid.major=element_blank(),#element_line(color='#dddddd',linetype="dashed",size=0.5),
  panel.background=element_rect(fill="white",colour='white'))+
  theme(legend.position='none',axis.text=element_blank())+
  guides(colour=guide_legend(override.aes=list(size=3)))+labs(tag="A")


#create barplot of morph frequencies
morph<-read.table("ordered_sample_info.txt",h=T)
morph$Cw_pheno<-factor(morph$Cw_pheno)
morph$Fw_pheno<-factor(morph$Fw_pheno)
morph$Sw_pheno<-factor(morph$Sw_pheno)
levels(morph$Cw_pheno)<-c("Straight","Curly")
levels(morph$Fw_pheno)<-c("Normal","Flat")
levels(morph$Sw_pheno)<-c("Long","Small")
summary(factor(morph$Population))

#calculate Cw proportion across sites
morph$Cw<-"N"
morph[grep("Cw",morph$Pheno),]$Cw<-"Y"
cw.table<-data.frame(cbind(table(morph$Population,morph$Cw)))
cw.table$prop.cw<-cw.table$Y/(cw.table$N+cw.table$Y)

# #calc singing proportion
morph$singing<-"N"
morph[grep("wtNw",morph$Pheno),]$singing<-"Y"
singing.table<-data.frame(cbind(table(morph$Population,morph$singing)))
singing.table$prop.cw<-singing.table$Y/(singing.table$N+singing.table$Y)
singing.table$Population<-rownames(singing.table)

singing.table.region<-data.frame(cbind(table(morph$Pop1_region,morph$singing)))
singing.table.region$prop.singing<-singing.table.region$Y/(singing.table.region$N+singing.table.region$Y)
singing.table.region$Population<-rownames(singing.table.region)

cw.table$site<-c("Hawaii.CL","Hawaii.UH","Kauai.AS","Kauai.CG","Kauai.PK","Kauai.VL","Kauai.HC","Oahu.BYU","Oahu.AC","Oahu.CC","Oahu.KP")

levels(morph$Pheno)
unique(morph$Population)
unique(morph$Phenotype)
morph[morph$Sw_pheno=="Small",]$Cw_pheno<-NA
morph$Phenotype<-factor(paste(morph$Sw_pheno,morph$Cw_pheno,morph$Fw_pheno,sep="_"),
                        levels=c("Long_Curly_Flat","Long_Curly_Normal","Long_Straight_Flat",
                                 "Small_NA_Normal","Long_Straight_Normal"))
morph[morph$Population=="Oahu_BYU",]$Population="Oahu.BYU"
morph[morph$Population=="Hilo.CL",]$Population="Hawaii.CL"
morph[morph$Population=="Hilo.UH",]$Population="Hawaii.UH"
morph$Population<-factor(morph$Population,levels=c("Kauai.WC","Kauai.CG","Kauai.AS","Kauai.VL","Kauai.PK",
                                                   "Oahu.CC","Oahu.BYU","Oahu.AC","Oahu.KP",
                                                   "Hawaii.CL","Hawaii.UH"))
#make long straight nw yello
g.freqs<-ggplot(morph,aes(x=Population,fill=Phenotype))+science_theme+scale_fill_viridis(discrete=TRUE,option='E')+
  geom_bar(position='fill',alpha=1)+theme(panel.grid=element_blank(),axis.text.x=element_text(angle=90,hjust = 1,vjust = 0.5))+
  theme(legend.position='right',legend.title=element_blank(),axis.title=element_blank(),
        legend.key.size = unit(0.45, "cm"))+
  guides(shape = guide_legend(override.aes = list(size = 0.5)))+
  guides(color = guide_legend(override.aes = list(size = 0.5)))+
  labs(tag="B")


library(gridExtra)
ggsave('patch_map_morph.svg',dpi=1200,width=7.25,height=3,units='in',plot=
         grid.arrange(Map_Populations_crickets,
                      g.freqs+theme(plot.margin = margin(75, 5.5, 5.5, 5.5, "pt")),widths=c(1.25,1)))

