library(ggbeeswarm)
library(ggplot2)
library(viridis)
library(lme4)
library(car)
library(rptR)
library(dplyr)

flies<-read.csv('flytrapping.csv',h=T)
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
g.time<-ggplot(flies,aes(x=time_postsunset,y=flies_caught))+
  geom_rect(fill='#fcf3f2',xmin=0,xmax=120,ymin=0,ymax=1.5)+coord_cartesian(ylim=c(0,1.5))+
  theme_bw()+theme(panel.grid=element_blank(),legend.position='none')+
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
rep1<-rpt(formula=(log(meanfly.yn+1)~(1|site)),grname="site",data=data.frame(flies.sum.night),datatype = "Gaussian")
qqnorm(resid(rep1$mod))#resids look reasonably normally distributed
qqline(resid(rep1$mod))
shapiro.test(resid(rep1$mod))#as above
plot(rep1$mod)#could be heteroscedasticity
leveneTest(log(meanfly.yn+1) ~ factor(site), data = data.frame(flies.sum.night))#P = 0.10. Probably borderline 
plot(rep1)

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
g.silent<-ggplot(flies.sum.night,aes(x=site,y=log(fw.rel.fitness),colour=island,group=date))+
  geom_hline(yintercept=0,linetype='dashed')+
  theme_bw()+theme(panel.grid=element_blank(),axis.text.x=element_text(angle=45,hjust=1,vjust=1),
                   axis.title.x=element_blank(),
                   strip.background=element_rect(fill='#eeeeee'))+
  facet_grid(.~island,scales='free',space='free')+geom_quasirandom(alpha=0.75,width=0.1)+
  scale_y_continuous(breaks=seq(-1,2,0.2))+
  ylab("log rel. fitness of silent males")+#scale_y_continuous(transform = 'log2')+annotation_logticks(sides = 'l')+
  stat_summary(aes(group=site,fill=island),colour='#555555',shape=21,linewidth=0.5,
               fun.data=mean_se,
               fun.args = list(mult = 2))+
  #geom_boxplot(aes(group=site),alpha=0)+
  scale_colour_manual(values=c("#00c08b","#c77cff","#f8766d"))+
  scale_fill_manual(values=c("#00c08b","#c77cff","#f8766d"))

#


##for outputting selection file for LFMM
flies.sum<-flies %>% filter(playback=='song' & time_postsunset<120) %>% 
  group_by(site,island) %>% 
  dplyr::summarise(meanflies=mean(flies_caught),meanfly.yn=mean(fly.yn))

Nten = flies.sum$meanflies
Nfly = Nten*0.76

#probability of  being infected per night is...
#(only 61% of attacks end in infestation)
Pinf = 1 - (1 - 0.61)^(Nfly)
flies.sum$Pinf = 1 - (1 - 0.61)^(Nfly)

selection<-read.csv("~/Documents/StA/popgen/LFMM/selection.csv")
unique(selection$pop)
unique(flies.sum$site)
flies.sum$pop<-NA
flies.sum[flies.sum$site=="Hawaii.CL",]$pop<-"Hilo_Church_Lawn"
flies.sum[flies.sum$site=="Hawaii.UH",]$pop<-"Hilo_UH"
flies.sum[flies.sum$site=="Kauai.AS",]$pop<-"Kauai_Agricultural_Station"
flies.sum[flies.sum$site=="Kauai.CG",]$pop<-"Kauai_Common_Ground"
flies.sum[flies.sum$site=="Kauai.HC",]$pop<-"Kauai_Waioli_Church"
flies.sum[flies.sum$site=="Kauai.PK",]$pop<-"Kauai_Pono_Kai"
flies.sum[flies.sum$site=="Kauai.VL",]$pop<-"Kauai_Vacant_Lot"
flies.sum[flies.sum$site=="Oahu.AC",]$pop<-"Oahu_Astronomy_Center"
flies.sum[flies.sum$site=="Oahu.BYU",]$pop<-"Oahu_BYU"
flies.sum[flies.sum$site=="Oahu.CC",]$pop<-"Oahu_Community_Center"
flies.sum[flies.sum$site=="Oahu.KP",]$pop<-"Oahu_Kamilo"

selection.risk<-inner_join(selection,flies.sum[,c("Pinf","pop")])
selection.risk<-selection.risk[,-3]
write.csv(selection.risk,file="~/Documents/StA/popgen/LFMM/selection_risk.csv",quote=FALSE,col.names = FALSE,row.names=FALSE)


#### investigate influence of rain and temperature
#load coordinates
sites<-read.csv("~/Documents/StA/popgen/map/cw_mapdata_silentvars.csv")
sites[sites$name1=="Kauai.WC",]$name1="Kauai.HC"

#read fly data
flies<-read.csv('~/Documents/StA/popgen/flytrapping/flytrapping.csv',h=T)
flies<-flies[flies$rain %in% c("dry","light rain"),]

flies.clim<-merge(flies.sum.night,sites,by.x='site',by.y='name1')
unique(flies.clim$site)

#retrieve daily temp
library(nasapower)
flies.clim$date1<-paste0(substr(flies.clim$date,7,11),"-",substr(flies.clim$date,4,5),"-",substr(flies.clim$date,1,2))
get_temp <- function(latitude, longitude, date) {
  result <- get_power(
    community = "AG",
    lonlat = c(longitude, latitude),
    pars = c("T2M"),  # T2M = Temperature at 2 meters
    dates = date,
    temporal_api = "daily"
  )
  print(result$T2M[1])
  return(result$T2M[1])
}

#retrieve daily rain
get_rain <- function(latitude, longitude, date) {
  result <- get_power(
    community = "AG",
    lonlat = c(longitude, latitude),
    pars = "PRECTOTCORR", 
    dates = date,
    temporal_api = "daily"
  )
  print(result$PRECTOTCORR[1])
  return(result$PRECTOTCORR[1])
}

#apply functions
flies.clim$temperature <- mapply(get_temp, flies.clim$latitude, flies.clim$longitude, flies.clim$date1)
flies.clim$rain <- mapply(get_rain, flies.clim$latitude, flies.clim$longitude, flies.clim$date1)
unique(flies.clim$site)

#glm.temp<-lmer(log(meanflies+1)~temperature+rain+(1|site),data=flies.clim)
glm.temp<-lm(log(meanflies+1)~temperature+rain+island.x/site,data=flies.clim)
qqnorm(resid(glm.temp))
qqline(resid(glm.temp))
car::Anova(glm.temp)
summary(glm.temp)

ggplot(flies.clim,aes(x=log(rain+1),y=log(meanflies.clim+1),colour=site))+
  theme_bw()+geom_point()+geom_smooth(method='lm',se=FALSE)+
  ylab("Daily mean flies.clim captured")+
  xlab("Log2 daily precipitation")#+facet_grid(.~temperature>median(temperature))


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
#   [1] nasapower_4.2.5   dplyr_1.1.4       rptR_0.9.22       car_3.1-3         carData_3.0-5    
# [6] lme4_1.1-36       Matrix_1.7-2      viridis_0.6.5     viridisLite_0.4.2 ggbeeswarm_0.7.2 
# [11] ggplot2_4.0.2    
# 
# loaded via a namespace (and not attached):
#   [1] gtable_0.3.6       compiler_4.4.2     tidyselect_1.2.1   Rcpp_1.1.0         gridExtra_2.3     
# [6] splines_4.4.2      scales_1.4.0       boot_1.3-31        lattice_0.22-6     R6_2.5.1          
# [11] generics_0.1.3     Formula_1.2-5      rbibutils_2.3      MASS_7.3-64        tibble_3.2.1      
# [16] nloptr_2.1.1       minqa_1.2.8        pillar_1.10.1      RColorBrewer_1.1-3 rlang_1.1.5       
# [21] S7_0.2.1           cli_3.6.3          withr_3.0.2        magrittr_2.0.3     Rdpack_2.6.2      
# [26] grid_4.4.2         rstudioapi_0.17.1  beeswarm_0.4.0     lifecycle_1.0.4    nlme_3.1-167      
# [31] vipor_0.4.7        reformulas_0.4.0   vctrs_0.6.5        glue_1.8.0         farver_2.1.2      
# [36] abind_1.4-8        tools_4.4.2        pkgconfig_2.0.3   

