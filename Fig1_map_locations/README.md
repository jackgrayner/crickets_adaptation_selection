# Fig1: map, locations

### cw_mapdata_silentvars.csv

Cricket population data: 
- island: Hawaiian island
- name2: location name
- name1: location name in final format used in paper
- singing_remain: whether singing males remain in population
- Phenotypes1: mutant wing phenotypes present in population
- Phenotypes: mutant wing phenotypes present in population (different formatting)
- longitude: longitude
- latitute: latitude

### ordered_sample_info.txt

Cricket genome sample sample data:
- Male_ID: assigned ID for animal
- SID: Sample ID used in read name
- Pheno: male assigned phenotype
- Population: Population in which male collected
- Pop1_region: Region in which male collected (if populations nearby)
- Cw_pheno: whether male expressed Cw (1=no, 2=yes)
- Fw_pheno: whether male expressed Fw (1=no, 2=yes)
- Sw_pheno: whether male expressed Sw (1=no, 2=yes)

### Plot_fig1.R

R script used to plot figure 1.

R version 4.4.2 (2024-10-31)
Platform: aarch64-apple-darwin20
Running under: macOS Sonoma 14.3

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: America/New_York
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] ggrepel_0.9.6           ape_5.8-1               wesanderson_0.3.7       viridis_0.6.5          
[5] viridisLite_0.4.2       ggspatial_1.1.9         rnaturalearthdata_1.0.0 rnaturalearth_1.0.1    
[9] ggplot2_4.0.2          

loaded via a namespace (and not attached):
 [1] gtable_0.3.6       jsonlite_1.8.9     dplyr_1.1.4        compiler_4.4.2     tidyselect_1.2.1  
 [6] Rcpp_1.1.0         parallel_4.4.2     gridExtra_2.3      scales_1.4.0       lattice_0.22-6    
[11] R6_2.5.1           generics_0.1.3     classInt_0.4-11    sf_1.0-19          tibble_3.2.1      
[16] units_0.8-5        DBI_1.2.3          pillar_1.10.1      RColorBrewer_1.1-3 rlang_1.1.5       
[21] terra_1.8-21       S7_0.2.1           cli_3.6.3          withr_3.0.2        magrittr_2.0.3    
[26] class_7.3-23       digest_0.6.37      grid_4.4.2         rstudioapi_0.17.1  nlme_3.1-167      
[31] lifecycle_1.0.4    vctrs_0.6.5        KernSmooth_2.23-26 proxy_0.4-27       glue_1.8.0        
[36] farver_2.1.2       codetools_0.2-20   e1071_1.7-16       httr_1.4.7         tools_4.4.2       
[41] pkgconfig_2.0.3   