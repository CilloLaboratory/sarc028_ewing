Introduction
------------

Here, we will use the previously defined cell types from scRNAseq
analysis to annotate the CosMx cell types.

Load packages
-------------

    .libPaths("/ix1/acillo/arc85/packages/Rlibs_sarc028")

    suppressMessages({
        library(Seurat)
        library(tidyverse)
        library(patchwork)
    })

Load datasets
-------------

    ser_scrnaseq <- readRDS("../../03_output/Ewing_scRNAseq_reference.rds")
    ser_cos <- readRDS("../../03_output/SARC028_CosMx_object.rds")

Combine scRNAseq and CosMx
--------------------------

    ## Select features present in both scRNAseq and Nanostring
    DefaultAssay(ser_scrnaseq) <- "RNA"

    sc_raw <- GetAssayData(ser_scrnaseq)

    DefaultAssay(ser_cos) <- "Nanostring"

    nano_raw <- GetAssayData(ser_cos)

    features_use <- intersect(rownames(sc_raw),rownames(nano_raw))

    sc_raw_fil <- sc_raw[features_use,]

    nano_raw_fil <- nano_raw[features_use,]

    ## Create new Seurat objects with only shared features
    sc_ser <- CreateSeuratObject(sc_raw_fil)

    nano_ser <- CreateSeuratObject(nano_raw_fil)

    sc_ser@meta.data <- ser_scrnaseq@meta.data

    sc_ser@meta.data$data_type <- "scrnaseq"

    nano_ser@meta.data <- ser_cos@meta.data

    nano_ser@meta.data$data_type <- "nanostring"

Integrate scRNAseq and CosMx data
---------------------------------

    ser_list <- list(sc_ser,nano_ser)

    # normalize and identify variable features for each dataset independently
    ser_list <- lapply(X = ser_list, FUN = function(x) {
        x <- NormalizeData(x)
        x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 500)
    })

    # select features that are repeatedly variable across datasets for integration
    features <- SelectIntegrationFeatures(object.list = ser_list)

    anchors <- FindIntegrationAnchors(object.list = ser_list, anchor.features = features)

    ## Scaling features for provided objects

    ## Finding all pairwise anchors

    ## Running CCA

    ## Merging objects

    ## Finding neighborhoods

    ## Finding anchors

    ##  Found 26356 anchors

    ## Filtering anchors

    ##  Retained 1464 anchors

    # this command creates an 'integrated' data assay
    ser_comb <- IntegrateData(anchorset = anchors)

    ## Merging dataset 1 into 2

    ## Extracting anchors for merged samples

    ## Finding integration vectors

    ## Finding integration vector weights

    ## Integrating data

    # We will perform downstream analysis on the corrected data
    # Original unmodified data still resides in the 'RNA' assay
    DefaultAssay(ser_comb) <- "integrated"

    # Run the standard workflow for visualization and clustering
    ser_comb <- ScaleData(ser_comb, verbose = FALSE)
    ser_comb <- RunPCA(ser_comb, npcs = 30, verbose = FALSE)

    ElbowPlot(ser_comb,ndims=50)

    ## Warning in ElbowPlot(ser_comb, ndims = 50): The object only has information for
    ## 30 reductions

![](03_integrated_cosmx_scrnaseq-exported_files/figure-markdown_strict/unnamed-chunk-4-1.png)

    ser_comb <- RunUMAP(ser_comb, reduction = "pca", dims = 1:20)

    ## Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
    ## To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
    ## This message will be shown once per session

    ## 10:19:00 UMAP embedding parameters a = 0.9922 b = 1.112

    ## 10:19:00 Read 63811 rows and found 20 numeric columns

    ## 10:19:00 Using Annoy for neighbor search, n_neighbors = 30

    ## 10:19:01 Building Annoy index with metric = cosine, n_trees = 50

    ## 0%   10   20   30   40   50   60   70   80   90   100%

    ## [----|----|----|----|----|----|----|----|----|----|

    ## **************************************************|
    ## 10:19:11 Writing NN index file to temp file /tmp/RtmpjvQNBv/filee6ed68750ee5
    ## 10:19:11 Searching Annoy index using 1 thread, search_k = 3000
    ## 10:19:52 Annoy recall = 99.69%
    ## 10:19:53 Commencing smooth kNN distance calibration using 1 thread with target n_neighbors = 30
    ## 10:19:55 295 smooth knn distance failures
    ## 10:19:57 Initializing from normalized Laplacian + noise (using irlba)
    ## 10:20:00 Commencing optimization for 200 epochs, with 3035836 positive edges
    ## 10:22:11 Optimization finished

    ser_comb <- FindNeighbors(ser_comb, reduction = "pca", dims = 1:20)

    ## Computing nearest neighbor graph
    ## Computing SNN

    ser_comb <- FindClusters(ser_comb, resolution = c(0.3,0.5,0.7,1))

    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 63811
    ## Number of edges: 1641797
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.9310
    ## Number of communities: 13
    ## Elapsed time: 19 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 63811
    ## Number of edges: 1641797
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.9107
    ## Number of communities: 18
    ## Elapsed time: 16 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 63811
    ## Number of edges: 1641797
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8943
    ## Number of communities: 19
    ## Elapsed time: 16 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 63811
    ## Number of edges: 1641797
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8715
    ## Number of communities: 26
    ## Elapsed time: 15 seconds

Check results
-------------

    DimPlot(ser_comb,group.by="data_type")

![](03_integrated_cosmx_scrnaseq-exported_files/figure-markdown_strict/unnamed-chunk-5-1.png)

    p1 <- DimPlot(ser_comb,group.by="scrnaseq_crude_cell_types",label=T,repel=T)
    p2 <- DimPlot(ser_comb,group.by="integrated_snn_res.0.5",label=T)

    p1 

![](03_integrated_cosmx_scrnaseq-exported_files/figure-markdown_strict/unnamed-chunk-5-2.png)

    p2 

![](03_integrated_cosmx_scrnaseq-exported_files/figure-markdown_strict/unnamed-chunk-5-3.png)

    Idents(ser_comb) <- "integrated_snn_res.0.5"

    ser_comb@meta.data %>%
        select(integrated_snn_res.0.5,nFeature_Nanostring) %>%
        filter(!is.na(nFeature_Nanostring)) %>%
        group_by(integrated_snn_res.0.5) %>%
        summarize(mean_features=mean(nFeature_Nanostring))

    ## # A tibble: 18 × 2
    ##    integrated_snn_res.0.5 mean_features
    ##    <fct>                          <dbl>
    ##  1 0                              42.7 
    ##  2 1                             120.  
    ##  3 2                              11.8 
    ##  4 3                             109.  
    ##  5 4                              59.7 
    ##  6 5                             122.  
    ##  7 6                              61.4 
    ##  8 7                              47.4 
    ##  9 8                              65.2 
    ## 10 9                              59.4 
    ## 11 10                             73.2 
    ## 12 11                             76.8 
    ## 13 12                            121.  
    ## 14 13                             84   
    ## 15 14                            107.  
    ## 16 15                             63   
    ## 17 16                              2.04
    ## 18 17                            173.

    ## Remove cluster 16 - essentially no features per cell
    ser_comb <- ser_comb[,!ser_comb@meta.data$integrated_snn_res.0.5=="16"]

    ser_comb@meta.data$integrated_snn_res.0.5 <- droplevels(ser_comb@meta.data$integrated_snn_res.0.5)

    table(ser_comb@meta.data$integrated_snn_res.0.5)

    ## 
    ##     0     1     2     3     4     5     6     7     8     9    10    11    12 
    ## 12992 12940  7758  4549  3846  3512  3357  2798  2479  2320  2268  1417  1230 
    ##    13    14    15    17 
    ##  1045   435   324   258

Transfer cell type annotations
------------------------------

    ser_comb@meta.data %>%
        select(integrated_snn_res.0.5,scrnaseq_crude_cell_types) %>%
        table() %>%
        pheatmap::pheatmap(.,scale="column")

![](03_integrated_cosmx_scrnaseq-exported_files/figure-markdown_strict/unnamed-chunk-6-1.png)

    meta <- ser_comb@meta.data %>%
        mutate(initial_cell_annot=recode(integrated_snn_res.0.5,
                                       `0` = "Tumor_cells_c0",
                                       `1` = "Tumor_cells_c1",
                                       `2` = "Tumor_cells_c2",
                                       `3` = "Tumor_cells_c3",
                                       `4` = "Osteoclasts",
                                       `5` = "Tumor_cells_c4",
                                       `6` = "T_and_NK_cells",
                                       `7` = "Tumor_cells_c5",
                                       `8` = "Tumor_cells_c6",
                                       `9` = "CD14posCD16pos_macro",
                                       `10` = "Tumor_cells_c7",
                                       `11` = "Monocytes",
                                       `12` = "CD1C_DCs",
                                       `13` = "Plasmablasts/pDCs",
                                       `14` = "Fibroblasts",
                                       `15` = "B_cells",
                                       `17` = "C1QCpos_macrophage")
               )

    table(meta$initial_cell_annot,useNA="always")

    ## 
    ##       Tumor_cells_c0       Tumor_cells_c1       Tumor_cells_c2 
    ##                12992                12940                 7758 
    ##       Tumor_cells_c3          Osteoclasts       Tumor_cells_c4 
    ##                 4549                 3846                 3512 
    ##       T_and_NK_cells       Tumor_cells_c5       Tumor_cells_c6 
    ##                 3357                 2798                 2479 
    ## CD14posCD16pos_macro       Tumor_cells_c7            Monocytes 
    ##                 2320                 2268                 1417 
    ##             CD1C_DCs    Plasmablasts/pDCs          Fibroblasts 
    ##                 1230                 1045                  435 
    ##              B_cells   C1QCpos_macrophage                 <NA> 
    ##                  324                  258                    0

    ser_comb@meta.data$initial_cell_annot <- meta$initial_cell_annot

    DimPlot(ser_comb,group.by="initial_cell_annot")

![](03_integrated_cosmx_scrnaseq-exported_files/figure-markdown_strict/unnamed-chunk-6-2.png)

Save results and session info
-----------------------------

    # saveRDS(ser_comb,file="../../03_output/SARC028_CosMx_and_scRNAseq_object.rds")

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
    ## [1] en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] patchwork_1.1.2    lubridate_1.9.2    forcats_1.0.0      stringr_1.5.0     
    ##  [5] dplyr_1.1.0        purrr_1.0.1        readr_2.1.4        tidyr_1.3.0       
    ##  [9] tibble_3.1.8       ggplot2_3.4.1      tidyverse_2.0.0    SeuratObject_4.1.3
    ## [13] Seurat_4.3.0.1     rmarkdown_2.25    
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] Rtsne_0.16             colorspace_2.0-3       deldir_1.0-6          
    ##   [4] ellipsis_0.3.2         ggridges_0.5.4         spatstat.data_3.0-0   
    ##   [7] farver_2.1.1           leiden_0.4.3           listenv_0.8.0         
    ##  [10] ggrepel_0.9.2          fansi_1.0.3            codetools_0.2-18      
    ##  [13] splines_4.2.0          knitr_1.39             polyclip_1.10-4       
    ##  [16] jsonlite_1.8.4         ica_1.0-3              cluster_2.1.3         
    ##  [19] png_0.1-7              pheatmap_1.0.12        uwot_0.1.14           
    ##  [22] shiny_1.7.3            sctransform_0.3.5      spatstat.sparse_3.0-0 
    ##  [25] compiler_4.2.0         httr_1.4.4             Matrix_1.5-3          
    ##  [28] fastmap_1.1.0          lazyeval_0.2.2         cli_3.6.0             
    ##  [31] later_1.3.0            htmltools_0.5.6.1      tools_4.2.0           
    ##  [34] igraph_1.3.5           gtable_0.3.1           glue_1.6.2            
    ##  [37] RANN_2.6.1             reshape2_1.4.4         Rcpp_1.0.9            
    ##  [40] scattermore_0.8        vctrs_0.5.2            spatstat.explore_3.0-5
    ##  [43] nlme_3.1-157           progressr_0.11.0       lmtest_0.9-40         
    ##  [46] spatstat.random_3.0-1  xfun_0.40              globals_0.16.2        
    ##  [49] timechange_0.2.0       mime_0.12              miniUI_0.1.1.1        
    ##  [52] lifecycle_1.0.3        irlba_2.3.5.1          goftest_1.2-3         
    ##  [55] future_1.29.0          MASS_7.3-57            zoo_1.8-11            
    ##  [58] scales_1.3.0           hms_1.1.2              promises_1.2.0.1      
    ##  [61] spatstat.utils_3.0-1   parallel_4.2.0         RColorBrewer_1.1-3    
    ##  [64] yaml_2.3.6             reticulate_1.26        pbapply_1.6-0         
    ##  [67] gridExtra_2.3          stringi_1.7.8          highr_0.9             
    ##  [70] rlang_1.1.1            pkgconfig_2.0.3        matrixStats_0.63.0    
    ##  [73] evaluate_0.15          lattice_0.20-45        ROCR_1.0-11           
    ##  [76] tensor_1.5             labeling_0.4.2         htmlwidgets_1.5.4     
    ##  [79] cowplot_1.1.1          tidyselect_1.2.0       parallelly_1.32.1     
    ##  [82] RcppAnnoy_0.0.20       plyr_1.8.8             magrittr_2.0.3        
    ##  [85] R6_2.5.1               generics_0.1.3         pillar_1.8.1          
    ##  [88] withr_2.5.0            fitdistrplus_1.1-8     survival_3.5-3        
    ##  [91] abind_1.4-5            sp_1.5-1               future.apply_1.10.0   
    ##  [94] KernSmooth_2.23-20     utf8_1.2.2             spatstat.geom_3.0-3   
    ##  [97] plotly_4.10.1          tzdb_0.3.0             grid_4.2.0            
    ## [100] data.table_1.14.6      digest_0.6.30          xtable_1.8-4          
    ## [103] httpuv_1.6.6           munsell_0.5.0          viridisLite_0.4.1
