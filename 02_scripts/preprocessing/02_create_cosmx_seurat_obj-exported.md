Introduction
------------

Here, we will utilize Seurat to create an object containing the SARC028
CosMx data.

Load packages
-------------

    .libPaths("/ix1/acillo/arc85/packages/Rlibs_sarc028")

    suppressMessages({
        library(Seurat)
        library(tidyverse)
        library(patchwork)
        library(here)
    })

Read in CosMx Data
------------------

    ser_cos <- LoadNanostring(data.dir = "../../01_input/SARC028_CosMx",fov="1")

    ## Warning: Invalid name supplied, making object name syntactically valid. New
    ## object name is X1; see ?make.names for more details on syntax validity

Load compiled metadata and clinical data
----------------------------------------

    metadata_to_add <- readRDS("../../01_input/SARC028_CosMx_deident_metadata.rds")

Add metadata to Seurat obj
--------------------------

    identical(metadata_to_add$cb,colnames(ser_cos))

    ## [1] TRUE

    ser_cos@meta.data$fov_selected <- metadata_to_add$fov_selected
    ser_cos@meta.data$deident_id <- metadata_to_add$deident_id
    ser_cos@meta.data$timepoint <- metadata_to_add$timepoint
    ser_cos@meta.data$tumor_perc_change <- metadata_to_add$tumor_perc_change
    ser_cos@meta.data$log_tumor_perc_change <- metadata_to_add$log_tumor_perc_change
    ser_cos@meta.data$best_resp <- metadata_to_add$best_resp

Normalize Data with SCT
-----------------------

    ser_cos <- SCTransform(ser_cos, assay = "Nanostring", clip.range = c(-10, 10))

    ## Calculating cell attributes from input UMI matrix: log_umi

    ## Variance stabilizing transformation of count matrix of size 978 by 50631

    ## Model formula is y ~ log_umi

    ## Get Negative Binomial regression parameters per gene

    ## Using 978 genes, 5000 cells

    ##   |                                                                              |                                                                      |   0%  |                                                                              |===================================                                   |  50%  |                                                                              |======================================================================| 100%

    ## Second step: Get residuals using fitted parameters for 978 genes

    ##   |                                                                              |                                                                      |   0%  |                                                                              |===================================                                   |  50%  |                                                                              |======================================================================| 100%

    ## Computing corrected count matrix for 978 genes

    ##   |                                                                              |                                                                      |   0%  |                                                                              |===================================                                   |  50%  |                                                                              |======================================================================| 100%

    ## Calculating gene attributes

    ## Wall clock passed: Time difference of 48.81713 secs

    ## Determine variable features

    ## Place corrected count matrix in counts slot

    ## Centering data matrix

    ## Set default assay to SCT

Save and session info
---------------------

    # saveRDS(ser_cos,file="../../03_output/SARC028_CosMx_object.rds")
    sessionInfo()

    ## R version 4.2.0 (2022-04-22)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Red Hat Enterprise Linux
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib64/libblas.so.3.4.2
    ## LAPACK: /usr/lib64/liblapack.so.3.4.2
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] here_1.0.1         patchwork_1.1.2    lubridate_1.9.2    forcats_1.0.0     
    ##  [5] stringr_1.4.1      dplyr_1.1.0        purrr_1.0.1        readr_2.1.4       
    ##  [9] tidyr_1.3.0        tibble_3.1.8       ggplot2_3.4.1      tidyverse_2.0.0   
    ## [13] SeuratObject_4.1.3 Seurat_4.3.0.1    
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] Rtsne_0.16             colorspace_2.0-3       deldir_1.0-6          
    ##   [4] ellipsis_0.3.2         ggridges_0.5.4         rprojroot_2.0.3       
    ##   [7] spatstat.data_3.0-0    leiden_0.4.3           listenv_0.8.0         
    ##  [10] ggrepel_0.9.2          fansi_1.0.3            codetools_0.2-18      
    ##  [13] splines_4.2.0          knitr_1.39             polyclip_1.10-4       
    ##  [16] jsonlite_1.8.4         ica_1.0-3              cluster_2.1.3         
    ##  [19] png_0.1-7              uwot_0.1.14            shiny_1.7.3           
    ##  [22] sctransform_0.3.5      spatstat.sparse_3.0-0  compiler_4.2.0        
    ##  [25] httr_1.4.4             Matrix_1.5-3           fastmap_1.1.0         
    ##  [28] lazyeval_0.2.2         cli_3.4.1              later_1.3.0           
    ##  [31] htmltools_0.5.6.1      tools_4.2.0            igraph_1.3.5          
    ##  [34] gtable_0.3.1           glue_1.6.2             RANN_2.6.1            
    ##  [37] reshape2_1.4.4         Rcpp_1.0.9             scattermore_0.8       
    ##  [40] vctrs_0.5.2            spatstat.explore_3.0-5 nlme_3.1-157          
    ##  [43] progressr_0.11.0       lmtest_0.9-40          spatstat.random_3.0-1 
    ##  [46] xfun_0.40              globals_0.16.2         timechange_0.2.0      
    ##  [49] mime_0.12              miniUI_0.1.1.1         lifecycle_1.0.3       
    ##  [52] irlba_2.3.5.1          goftest_1.2-3          future_1.29.0         
    ##  [55] MASS_7.3-57            zoo_1.8-11             scales_1.3.0          
    ##  [58] hms_1.1.2              promises_1.2.0.1       spatstat.utils_3.0-1  
    ##  [61] parallel_4.2.0         RColorBrewer_1.1-3     yaml_2.3.6            
    ##  [64] reticulate_1.26        pbapply_1.6-0          gridExtra_2.3         
    ##  [67] stringi_1.7.8          rlang_1.1.1            pkgconfig_2.0.3       
    ##  [70] matrixStats_0.63.0     evaluate_0.15          lattice_0.20-45       
    ##  [73] ROCR_1.0-11            tensor_1.5             htmlwidgets_1.5.4     
    ##  [76] cowplot_1.1.1          tidyselect_1.2.0       parallelly_1.32.1     
    ##  [79] RcppAnnoy_0.0.20       plyr_1.8.8             magrittr_2.0.3        
    ##  [82] R6_2.5.1               generics_0.1.3         pillar_1.8.1          
    ##  [85] withr_2.5.0            fitdistrplus_1.1-8     survival_3.5-3        
    ##  [88] abind_1.4-5            sp_1.5-1               future.apply_1.10.0   
    ##  [91] KernSmooth_2.23-20     utf8_1.2.2             spatstat.geom_3.0-3   
    ##  [94] plotly_4.10.1          tzdb_0.3.0             rmarkdown_2.25        
    ##  [97] grid_4.2.0             data.table_1.14.6      digest_0.6.30         
    ## [100] xtable_1.8-4           httpuv_1.6.6           munsell_0.5.0         
    ## [103] viridisLite_0.4.1
