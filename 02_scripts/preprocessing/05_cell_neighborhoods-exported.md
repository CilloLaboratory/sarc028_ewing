Introduction
------------

Identify spatially localized cell neighborhoods based on cell type
annotations and spatial localization from the CosMx data.

Load packages
-------------

    .libPaths("/ihome/acillo/arc85/Rlibs_Mar_2023")

    suppressMessages({
        library(Seurat)
        library(tidyverse)
        library(patchwork)
        library(cluster)
        library(parallel)
        library(here)
    })

Load datasets
-------------

    ser_cos <- readRDS("../../03_output/SARC028_CosMx_object.rds")

    ser_comb <- readRDS("../../03_output/SARC028_CosMx_and_scRNAseq_object_refined_annotations.rds")

Filter cells and include cell annotations
-----------------------------------------

    ## Excluded cells filtered in combined analysis
    comb_meta <- ser_comb@meta.data %>%
        as_tibble(.,rownames="cb") %>%
        filter(data_type=="nanostring") %>%
        select(cb,final_cell_annot)

    dim(ser_cos)

    ## [1]   978 50631

    ser_cos <- ser_cos[,comb_meta$cb]

    dim(ser_cos)

    ## [1]   978 50348

    ## Include cell annotations
    meta <- ser_cos@meta.data %>%
        as_tibble(.,rownames="cb") %>%
        left_join(.,comb_meta,by="cb")

    table(meta$final_cell_annot)

    ## 
    ##              B_cells   C1QCpos_macrophage CD14posCD16pos_macro 
    ##                    1                    3                  601 
    ##             CD1C_DCs         CD4pos_Tconv          CD4pos_Treg 
    ##                  198                  392                    2 
    ##        CD8pos_T_cell          Fibroblasts            Monocytes 
    ##                   43                    5                   44 
    ##              NK_cell          Osteoclasts                 pDCs 
    ##                   28                 3589                    1 
    ##         Plasmablasts       Tumor_cells_c0       Tumor_cells_c1 
    ##                   15                12726                12885 
    ##       Tumor_cells_c2       Tumor_cells_c3       Tumor_cells_c4 
    ##                 7618                 4454                  893 
    ##       Tumor_cells_c5       Tumor_cells_c6       Tumor_cells_c7 
    ##                 2794                 1968                 2088

    ser_cos@meta.data$final_cell_annot <- meta$final_cell_annot

    ## Combine B lineage and filter rare cell types
    ser_cos@meta.data$final_cell_annot <- ifelse(ser_cos@meta.data$final_cell_annot %in% c("B_cells","Plasmablasts"),"B_cell_lineage",ser_cos@meta.data$final_cell_annot)

    cell_types_include <- table(ser_cos@meta.data$final_cell_annot) %>%
        data.frame() %>%
        filter(Freq>9) %>%
        pull(Var1)

    ser_cos <- ser_cos[,ser_cos@meta.data$final_cell_annot %in% cell_types_include]

    table(ser_cos@meta.data$final_cell_annot)

    ## 
    ##       B_cell_lineage CD14posCD16pos_macro             CD1C_DCs 
    ##                   16                  601                  198 
    ##         CD4pos_Tconv        CD8pos_T_cell            Monocytes 
    ##                  392                   43                   44 
    ##              NK_cell          Osteoclasts       Tumor_cells_c0 
    ##                   28                 3589                12726 
    ##       Tumor_cells_c1       Tumor_cells_c2       Tumor_cells_c3 
    ##                12885                 7618                 4454 
    ##       Tumor_cells_c4       Tumor_cells_c5       Tumor_cells_c6 
    ##                  893                 2794                 1968 
    ##       Tumor_cells_c7 
    ##                 2088

Function to find vector of neighbors
------------------------------------

Define a function to identify the cells with a specific radius of other
cells.

    cell_neighbor <- function(seurat_image,cell_types) {

        dat <- FetchData(seurat_image,vars="centroids")
        dat$cell_types <- cell_types

        dist_dat <- as.matrix(dist(dat[,1:2],method="euclidian"))

        neighbor_list <- vector("list",length=nrow(dat))

        for (i in 1:nrow(dat)) { 

            dat_sub <- dat
            dat_sub$cell_dist <- dist_dat[,i]

            dat_sub$include <- dat_sub$cell_dist<300

            neighbor_list[[i]] <- dat_sub[-i,] %>%
                filter(include==TRUE) %>%
                select(cell_types) %>%
                table() %>%
                prop.table() %>%
                as.vector()

        }

        neighbor_frame <- do.call(rbind,neighbor_list)
        colnames(neighbor_frame) <- levels(dat$cell_types)

        return(neighbor_frame)
        
    }

Cell neighborhoods across all images
------------------------------------

    ## Make sure cell types are a factor
    ser_cos@meta.data$final_cell_annot <- as.factor(ser_cos@meta.data$final_cell_annot)

    ## Split into individual FOVs
    ser_cos_split <- SplitObject(ser_cos,split="fov_selected")

    ## Loop over FOVs to ID cell neighorhoods
    ## Use parallel processing
    set.seed(250306)
    hoods <- mclapply(ser_cos_split,mc.cores=4,mc.set.seed=TRUE,function(x) {
        cell_neighbor(x$X1,cell_types=x@meta.data$final_cell_annot)
    })

    ## Combine into dataframe
    hood_frame_list <- vector("list",length=length(hoods))

    for (i in 1:length(hoods)) {
        hood_frame_list[[i]] <- data.frame(hoods[[i]],fov=names(ser_cos_split)[[i]])
    }

    hoods_frame <- do.call(rbind,hood_frame_list)

Evaluate spatial clustering solutions
-------------------------------------

    ## Remove some cells that were not near other cells 
    rows_exclude <- apply(hoods_frame,1,function(x) any(is.na(x)))

    ## Keep track of cells that get excluded
    hoods_frame$cb <- colnames(ser_cos)

    hoods_frame <- hoods_frame[!rows_exclude,]

    ## Evaluate clustering solutions
    hoods_res <- vector("list",length=7)
    hoods_dist <- dist(hoods_frame[,1:16])

    for (i in 1:length(hoods_res)) {
        clust_res <- data.frame(hoods_frame,cluster=kmeans(hoods_frame[,1:16],centers=3+i)$cluster)
        hoods_res[[i]] <- silhouette(clust_res$cluster,dist=hoods_dist)
        print(paste("Finished ",i," of 7",sep=""))
    }

    ## [1] "Finished 1 of 7"
    ## [1] "Finished 2 of 7"
    ## [1] "Finished 3 of 7"
    ## [1] "Finished 4 of 7"
    ## [1] "Finished 5 of 7"
    ## [1] "Finished 6 of 7"
    ## [1] "Finished 7 of 7"

    clustering_results <- data.frame(cluster_number=seq(4,10,1),avg_silhouette_width=sapply(hoods_res,function(x) summary(x)$avg.width))

    clustering_results %>%
        ggplot(.,aes(x=cluster_number,y=avg_silhouette_width)) +
        geom_col() +
        theme_bw()

![](05_cell_neighborhoods-exported_files/figure-markdown_strict/unnamed-chunk-6-1.png)

    ## choosing 10 clusters

Create spatial clustering with 10 clusters
------------------------------------------

    set.seed(250306)
    hoods_clust <- data.frame(hoods_frame,cluster=kmeans(hoods_frame[,1:16],centers=10)$cluster)

    mean_cell_freq <- hoods_clust %>%
        gather(cell_type,freq,-cluster,-fov,-cb) %>%
        group_by(cluster,cell_type) %>%
        summarize(mean_freq=mean(freq)) %>%
        ungroup() %>%
        spread(cell_type,mean_freq)

    ## `summarise()` has grouped output by 'cluster'. You can override using the
    ## `.groups` argument.

    mean_cell_freq <- mean_cell_freq[,2:ncol(mean_cell_freq)] %>%
        data.frame()

    pheatmap::pheatmap(mean_cell_freq,scale="column",color=colorRampPalette(colors = c("purple","black","yellow"))(1000))

![](05_cell_neighborhoods-exported_files/figure-markdown_strict/unnamed-chunk-7-1.png)

    ## Renaming clustering to match published version
    mapping <- data.frame(here_clust=seq(1,10,1),pub_clust=NA)

    mapping[3,"pub_clust"] <- 8
    mapping[9,"pub_clust"] <- 4
    mapping[2,"pub_clust"] <- 5
    mapping[5,"pub_clust"] <- 2
    mapping[1,"pub_clust"] <- 1
    mapping[10,"pub_clust"] <- 9
    mapping[4,"pub_clust"] <- 7
    mapping[7,"pub_clust"] <- 6
    mapping[8,"pub_clust"] <- 3
    mapping[6,"pub_clust"] <- 10

    rownames(mean_cell_freq) <- paste("Hood ",c(mapping$pub_clust),sep="")

    ## Plotting
    pheatmap::pheatmap(mean_cell_freq,scale="column",color=colorRampPalette(colors = c("purple","black","yellow"))(1000))

![](05_cell_neighborhoods-exported_files/figure-markdown_strict/unnamed-chunk-7-2.png)

    ## Rename clusters in hood_clust
    mapping

    ##    here_clust pub_clust
    ## 1           1         1
    ## 2           2         5
    ## 3           3         8
    ## 4           4         7
    ## 5           5         2
    ## 6           6        10
    ## 7           7         6
    ## 8           8         3
    ## 9           9         4
    ## 10         10         9

    hood_clust <- hoods_clust %>%
        mutate(cluster=recode(cluster,
                              `1`=1,
                              `2`=5,
                              `3`=8,
                              `4`=7,
                              `5`=2,
                              `6`=10,
                              `7`=6,
                              `8`=3,
                              `9`=4,
                              `10`=9)
               )

Add spatial clusters to Seurat object
-------------------------------------

    meta <- ser_cos@meta.data %>%
        as_tibble(.,rownames="cb")

    meta <- left_join(meta,hoods_clust,by="cb")

    ser_cos@meta.data$spatial_clusters <- meta$cluster

Find DEGs across spatial neighborhoods
--------------------------------------

    ser_cos <- ser_cos[,!is.na(ser_cos@meta.data$spatial_clusters)]

    Idents(ser_cos) <- "spatial_clusters"

    hoods_degs <- FindAllMarkers(ser_cos,logfc.threshold=0.1)

    ## Calculating cluster 3

    ## Calculating cluster 1

    ## Calculating cluster 5

    ## Calculating cluster 2

    ## Calculating cluster 10

    ## Calculating cluster 6

    ## Calculating cluster 8

    ## Calculating cluster 7

    ## Calculating cluster 4

    ## Calculating cluster 9

    degs_use <- hoods_degs %>%
        filter(p_val_adj<0.05) %>%
        filter(avg_log2FC>0) %>%
        select(gene) %>%
        distinct() %>%
        pull()

    dat_use <- GetAssayData(ser_cos,slot="data",assay="SCT") %>%
        as.matrix()

    ser_cos_genes <- cbind(ser_cos@meta.data,t(dat_use))

Save data and session info
--------------------------

    # saveRDS(hoods_clust,file="../../03_output/SARC028_spatial_clustering_results.rds")
    # saveRDS(hoods_degs,file="../../03_output/SARC028_CosMx_spatial_hood_DEGs.rds")
    # saveRDS(ser_cos_genes,file="../../03_output/SARC028_CosMx_expression_metadata.rds")

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
    ## [1] parallel  stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] here_1.0.1         cluster_2.1.3      patchwork_1.1.2    lubridate_1.9.2   
    ##  [5] forcats_1.0.0      stringr_1.5.0      dplyr_1.1.0        purrr_1.0.1       
    ##  [9] readr_2.1.4        tidyr_1.3.0        tibble_3.1.8       ggplot2_3.4.1     
    ## [13] tidyverse_2.0.0    SeuratObject_4.1.3 Seurat_4.3.0.1     rmarkdown_2.25    
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] Rtsne_0.16             colorspace_2.0-3       deldir_1.0-6          
    ##   [4] ellipsis_0.3.2         ggridges_0.5.4         rprojroot_2.0.3       
    ##   [7] spatstat.data_3.0-0    farver_2.1.1           leiden_0.4.3          
    ##  [10] listenv_0.8.0          ggrepel_0.9.2          fansi_1.0.3           
    ##  [13] codetools_0.2-18       splines_4.2.0          knitr_1.39            
    ##  [16] polyclip_1.10-4        jsonlite_1.8.4         ica_1.0-3             
    ##  [19] png_0.1-7              pheatmap_1.0.12        uwot_0.1.14           
    ##  [22] shiny_1.7.3            sctransform_0.3.5      spatstat.sparse_3.0-0 
    ##  [25] compiler_4.2.0         httr_1.4.4             Matrix_1.5-3          
    ##  [28] fastmap_1.1.0          lazyeval_0.2.2         limma_3.52.4          
    ##  [31] cli_3.6.0              later_1.3.0            htmltools_0.5.6.1     
    ##  [34] tools_4.2.0            igraph_1.3.5           gtable_0.3.1          
    ##  [37] glue_1.6.2             RANN_2.6.1             reshape2_1.4.4        
    ##  [40] Rcpp_1.0.9             scattermore_0.8        vctrs_0.5.2           
    ##  [43] spatstat.explore_3.0-5 nlme_3.1-157           progressr_0.11.0      
    ##  [46] lmtest_0.9-40          spatstat.random_3.0-1  xfun_0.40             
    ##  [49] globals_0.16.2         timechange_0.2.0       mime_0.12             
    ##  [52] miniUI_0.1.1.1         lifecycle_1.0.3        irlba_2.3.5.1         
    ##  [55] goftest_1.2-3          future_1.29.0          MASS_7.3-57           
    ##  [58] zoo_1.8-11             scales_1.3.0           hms_1.1.2             
    ##  [61] promises_1.2.0.1       spatstat.utils_3.0-1   RColorBrewer_1.1-3    
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
