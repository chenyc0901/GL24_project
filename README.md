# GL24_project
The analysis for DE genes upon GL24 treatment 


R version 4.2.2 (2022-10-31 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 22000)

attached base packages:
[1] grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] heatmap3_1.1.9         EnhancedVolcano_1.16.0 ggrepel_0.9.1          factoextra_1.0.7       gplots_3.1.3          
 [6] VennDiagram_1.7.3      futile.logger_1.4.3    corrplot_0.92          enrichplot_1.18.0      biomaRt_2.54.0        
[11] edgeR_3.40.0           limma_3.54.0           forcats_0.5.2          stringr_1.4.1          dplyr_1.0.10          
[16] purrr_0.3.5            readr_2.1.3            tidyr_1.2.1            tibble_3.1.8           ggplot2_3.4.0         
[21] tidyverse_1.3.2       

loaded via a namespace (and not attached):
  [1] readxl_1.4.1           shadowtext_0.1.2       backports_1.4.1        fastmatch_1.1-3        BiocFileCache_2.6.0   
  [6] plyr_1.8.7             igraph_1.3.5           lazyeval_0.2.2         splines_4.2.2          BiocParallel_1.32.0   
 [11] GenomeInfoDb_1.34.3    digest_0.6.30          yulab.utils_0.0.5      GOSemSim_2.24.0        viridis_0.6.2         
 [16] GO.db_3.16.0           fansi_1.0.3            magrittr_2.0.3         memoise_2.0.1          googlesheets4_1.0.1   
 [21] tzdb_0.3.0             fastcluster_1.2.3      Biostrings_2.66.0      graphlayouts_0.8.3     modelr_0.1.10         
 [26] prettyunits_1.1.1      colorspace_2.0-3       blob_1.2.3             rvest_1.0.3            rappdirs_0.3.3        
 [31] haven_2.5.1            crayon_1.5.2           RCurl_1.98-1.9         jsonlite_1.8.3         scatterpie_0.1.8      
 [36] ape_5.6-2              glue_1.6.2             polyclip_1.10-4        gtable_0.3.1           gargle_1.2.1          
 [41] zlibbioc_1.44.0        XVector_0.38.0         car_3.1-1              BiocGenerics_0.44.0    abind_1.4-5           
 [46] scales_1.2.1           DOSE_3.24.0            futile.options_1.0.1   DBI_1.1.3              rstatix_0.7.1         
 [51] Rcpp_1.0.9             viridisLite_0.4.1      progress_1.2.2         gridGraphics_0.5-1     tidytree_0.4.1        
 [56] bit_4.0.4              stats4_4.2.2           httr_1.4.4             fgsea_1.24.0           RColorBrewer_1.1-3    
 [61] ellipsis_0.3.2         pkgconfig_2.0.3        XML_3.99-0.12          farver_2.1.1           dbplyr_2.2.1          
 [66] locfit_1.5-9.6         utf8_1.2.2             labeling_0.4.2         ggplotify_0.1.0        tidyselect_1.2.0      
 [71] rlang_1.0.6            reshape2_1.4.4         AnnotationDbi_1.60.0   munsell_0.5.0          cellranger_1.1.0      
 [76] tools_4.2.2            cachem_1.0.6           cli_3.4.1              generics_0.1.3         RSQLite_2.2.18        
 [81] broom_1.0.1            fastmap_1.1.0          ggtree_3.6.2           bit64_4.0.5            fs_1.5.2              
 [86] tidygraph_1.2.2        caTools_1.18.2         KEGGREST_1.38.0        ggraph_2.1.0           nlme_3.1-160          
 [91] formatR_1.12           aplot_0.1.8            xml2_1.3.3             compiler_4.2.2         rstudioapi_0.14       
 [96] filelock_1.0.2         curl_4.3.3             png_0.1-7              ggsignif_0.6.4         reprex_2.0.2          
[101] treeio_1.22.0          tweenr_2.0.2           stringi_1.7.8          lattice_0.20-45        Matrix_1.5-1          
[106] vctrs_0.5.0            pillar_1.8.1           lifecycle_1.0.3        data.table_1.14.4      cowplot_1.1.1         
[111] bitops_1.0-7           patchwork_1.1.2        qvalue_2.30.0          R6_2.5.1               KernSmooth_2.23-20    
[116] gridExtra_2.3          IRanges_2.32.0         codetools_0.2-18       lambda.r_1.2.4         gtools_3.9.3          
[121] MASS_7.3-58.1          assertthat_0.2.1       withr_2.5.0            S4Vectors_0.36.0       GenomeInfoDbData_1.2.9
[126] parallel_4.2.2         hms_1.1.2              ggfun_0.0.8            HDO.db_0.99.1          carData_3.0-5         
[131] googledrive_2.0.0      ggpubr_0.5.0           ggforce_0.4.1          Biobase_2.58.0         lubridate_1.8.0 