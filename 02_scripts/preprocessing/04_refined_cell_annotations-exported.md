Introduction
------------

Here, we will perform subclustering to identify finer subsets of immune
populations within the CosMx data.

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

    ser_comb <- readRDS("../../03_output/SARC028_CosMx_and_scRNAseq_object.rds")

    DimPlot(ser_comb,group.by="initial_cell_annot",label=T)

![](04_refined_cell_annotations-exported_files/figure-markdown_strict/unnamed-chunk-2-1.png)

Isolate data for NK and T cell subset identification
----------------------------------------------------

    ser_comb <- ser_comb[,grepl("T_and_NK_cells",ser_comb@meta.data$initial_cell_annot)]

    table(ser_comb@meta.data$initial_cell_annot)

    ## 
    ##       Tumor_cells_c0       Tumor_cells_c1       Tumor_cells_c2 
    ##                    0                    0                    0 
    ##       Tumor_cells_c3          Osteoclasts       Tumor_cells_c4 
    ##                    0                    0                    0 
    ##       T_and_NK_cells       Tumor_cells_c5       Tumor_cells_c6 
    ##                 3357                    0                    0 
    ## CD14posCD16pos_macro       Tumor_cells_c7            Monocytes 
    ##                    0                    0                    0 
    ##             CD1C_DCs    Plasmablasts/pDCs          Fibroblasts 
    ##                    0                    0                    0 
    ##              B_cells   C1QCpos_macrophage 
    ##                    0                    0

    DimPlot(ser_comb,group.by="initial_cell_annot")

![](04_refined_cell_annotations-exported_files/figure-markdown_strict/unnamed-chunk-3-1.png)

Combine and integrate CosMx and scRNAseq for NK and T cells
-----------------------------------------------------------

    DefaultAssay(ser_comb) <- "RNA"

    ser_split <- SplitObject(ser_comb,split.by="data_type")

    ser_list <- ser_split

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

    ##  Found 2312 anchors

    ## Filtering anchors

    ##  Retained 1808 anchors

    # this command creates an 'integrated' data assay
    ser_comb <- IntegrateData(anchorset = anchors)

    ## Merging dataset 2 into 1

    ## Extracting anchors for merged samples

    ## Finding integration vectors

    ## Finding integration vector weights

    ## Integrating data

    # specify that we will perform downstream analysis on the corrected data note that the
    # original unmodified data still resides in the 'RNA' assay
    DefaultAssay(ser_comb) <- "integrated"

    # Run the standard workflow for visualization and clustering
    ser_comb <- ScaleData(ser_comb, verbose = FALSE)
    ser_comb <- RunPCA(ser_comb, npcs = 30, verbose = FALSE)

    ElbowPlot(ser_comb)

![](04_refined_cell_annotations-exported_files/figure-markdown_strict/unnamed-chunk-4-1.png)

    ser_comb <- RunUMAP(ser_comb, reduction = "pca", dims = 1:12)

    ## 10:51:19 UMAP embedding parameters a = 0.9922 b = 1.112

    ## 10:51:19 Read 3357 rows and found 12 numeric columns

    ## 10:51:19 Using Annoy for neighbor search, n_neighbors = 30

    ## 10:51:19 Building Annoy index with metric = cosine, n_trees = 50

    ## 0%   10   20   30   40   50   60   70   80   90   100%

    ## [----|----|----|----|----|----|----|----|----|----|

    ## **************************************************|
    ## 10:51:19 Writing NN index file to temp file /tmp/RtmpjvQNBv/filee6ed6579e8be
    ## 10:51:19 Searching Annoy index using 1 thread, search_k = 3000
    ## 10:51:21 Annoy recall = 100%
    ## 10:51:21 Commencing smooth kNN distance calibration using 1 thread with target n_neighbors = 30
    ## 10:51:22 Initializing from normalized Laplacian + noise (using irlba)
    ## 10:51:22 Commencing optimization for 500 epochs, with 137940 positive edges
    ## 10:51:38 Optimization finished

    ser_comb <- FindNeighbors(ser_comb, reduction = "pca", dims = 1:12)

    ## Computing nearest neighbor graph
    ## Computing SNN

    ser_comb <- FindClusters(ser_comb, resolution = c(0.3,0.5,0.7,1))

    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 3357
    ## Number of edges: 114094
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8793
    ## Number of communities: 8
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 3357
    ## Number of edges: 114094
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8423
    ## Number of communities: 10
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 3357
    ## Number of edges: 114094
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8157
    ## Number of communities: 11
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 3357
    ## Number of edges: 114094
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.7832
    ## Number of communities: 14
    ## Elapsed time: 0 seconds

Check results for NK and T cells
--------------------------------

    DimPlot(ser_comb,group.by="data_type")

![](04_refined_cell_annotations-exported_files/figure-markdown_strict/unnamed-chunk-5-1.png)

    p1 <- FeaturePlot(ser_comb,c("FOXP3","CD8A","FCGR3A","CD4"))

    p2 <- DimPlot(ser_comb,group.by="integrated_snn_res.0.5",label=T)

    p3 <- DimPlot(ser_comb,group.by="integrated_snn_res.1",label=T)

    p1 + p2

![](04_refined_cell_annotations-exported_files/figure-markdown_strict/unnamed-chunk-5-2.png)

    p3

![](04_refined_cell_annotations-exported_files/figure-markdown_strict/unnamed-chunk-5-3.png)

    Idents(ser_comb) <- "integrated_snn_res.1"

    cluster7 <- FindMarkers(ser_comb,ident.1="7")
    cluster8 <- FindMarkers(ser_comb,ident.1="8")

    cluster7 %>%
        filter(avg_log2FC>0) %>%
        filter(p_val_adj<0.05)

    ##                p_val avg_log2FC pct.1 pct.2     p_val_adj
    ## TYROBP 1.545677e-120  3.0013909 0.940 0.280 1.120616e-117
    ## FCGR3A 1.281933e-101  3.9946149 0.678 0.090  9.294017e-99
    ## GNLY   4.052346e-101  2.6966836 0.907 0.293  2.937951e-98
    ## FCER1G  5.769553e-99  3.2952369 0.809 0.229  4.182926e-96
    ## GZMB    3.826759e-98  2.7959165 0.880 0.293  2.774400e-95
    ## PRF1    2.495834e-88  2.6787082 0.869 0.338  1.809480e-85
    ## NKG7    2.802099e-82  1.7186122 0.973 0.606  2.031522e-79
    ## PLAC8   2.375975e-63  2.6153379 0.721 0.257  1.722582e-60
    ## KLRB1   8.320909e-62  1.7585493 0.918 0.438  6.032659e-59
    ## CST7    5.976088e-45  1.5059778 0.863 0.535  4.332664e-42
    ## IGFBP7  3.516378e-41  3.7926305 0.415 0.049  2.549374e-38
    ## GZMA    2.805709e-36  1.0904700 0.907 0.609  2.034139e-33
    ## CCL4    2.266546e-34  1.2735222 0.852 0.500  1.643246e-31
    ## CTSW    2.512189e-33  1.3050897 0.852 0.561  1.821337e-30
    ## CCL3    2.765724e-33  2.2090209 0.596 0.261  2.005150e-30
    ## FGR     3.980052e-28  2.8434760 0.399 0.124  2.885538e-25
    ## GZMH    4.457035e-19  1.3308462 0.612 0.355  3.231350e-16
    ## LAIR1   1.231850e-18  2.3087438 0.404 0.181  8.930913e-16
    ## ADGRG1  1.958289e-12  3.0437229 0.322 0.138  1.419759e-09
    ## XBP1    4.441842e-12  1.1757688 0.514 0.331  3.220335e-09
    ## CX3CR1  9.812900e-11  3.1195761 0.235 0.062  7.114352e-08
    ## TXK     9.271810e-10  1.8369821 0.350 0.192  6.722063e-07
    ## CD300A  1.536771e-09  1.9537480 0.355 0.195  1.114159e-06
    ## MALAT1  2.072726e-09  0.3041943 1.000 0.997  1.502726e-06
    ## NEAT1   2.844463e-06  0.5896483 0.678 0.597  2.062236e-03
    ## RAMP1   3.839565e-06  1.6993810 0.164 0.068  2.783684e-03
    ## XCL2    5.820023e-06  1.1920656 0.377 0.255  4.219516e-03
    ## HAVCR2  6.284751e-06  1.9849196 0.257 0.145  4.556444e-03
    ## NCR1    1.186151e-05  2.2642549 0.191 0.083  8.599594e-03

    cluster8 %>%
        filter(avg_log2FC>0) %>%
        filter(p_val_adj<0.05)

    ##                 p_val avg_log2FC pct.1 pct.2    p_val_adj
    ## HSPA1A   1.377668e-32  2.6826963 0.512 0.158 9.988091e-30
    ## NKG7     2.532042e-20  0.8178557 0.889 0.612 1.835731e-17
    ## HSPA1B   3.672442e-20  2.4382127 0.451 0.196 2.662520e-17
    ## MT2A     1.328951e-15  1.2109498 0.630 0.401 9.634896e-13
    ## GZMK     3.746555e-14  1.0811681 0.648 0.441 2.716252e-11
    ## CYTOR    4.418236e-14  1.9766545 0.426 0.219 3.203221e-11
    ## CD27     7.718904e-14  1.0757426 0.636 0.425 5.596205e-11
    ## CST7     1.037275e-13  0.7073038 0.765 0.542 7.520246e-11
    ## HLA-DRB1 2.575811e-13  0.8989496 0.691 0.461 1.867463e-10
    ## GZMA     4.975155e-10  0.5045323 0.778 0.618 3.606987e-07
    ## LGALS1   2.184283e-09  0.6967232 0.722 0.538 1.583605e-06
    ## CD8A     9.734837e-09  0.8238041 0.611 0.448 7.057757e-06
    ## IL2RG    4.858798e-08  0.2863562 0.784 0.628 3.522628e-05
    ## CTSW     8.728869e-08  0.5077513 0.722 0.570 6.328430e-05
    ## CCL4     1.343961e-07  0.6689457 0.660 0.512 9.743716e-05
    ## LAG3     3.352291e-07  1.6687030 0.370 0.240 2.430411e-04
    ## HLA-DPB1 3.717234e-07  0.6762718 0.654 0.495 2.694995e-04
    ## HLA-DQA1 5.658836e-06  1.5706990 0.327 0.214 4.102656e-03
    ## SELENOP  6.695755e-06  2.4326174 0.167 0.085 4.854422e-03
    ## ITM2A    8.334844e-06  0.5977685 0.617 0.488 6.042762e-03
    ## CXCR3    2.573387e-05  0.8041912 0.488 0.371 1.865706e-02
    ## VCAM1    3.206462e-05  2.9524116 0.123 0.050 2.324685e-02
    ## GZMH     6.717370e-05  0.8355907 0.475 0.364 4.870093e-02

Add refined metadata for NK and T cells
---------------------------------------

    meta <- ser_comb@meta.data %>%
        mutate(refined_lymph_cell_annot=recode(integrated_snn_res.1,
                                       `0` = "CD8pos_T_cell",
                                       `1` = "CD4pos_Tconv",
                                       `2` = "CD4pos_Tconv",
                                       `3` = "CD4pos_Tconv",
                                       `4` = "CD4pos_Tconv",
                                       `5` = "CD8pos_T_cell",
                                       `6` = "CD8pos_T_cell",
                                       `7` = "NK_cell",
                                       `8` = "CD8pos_T_cell",
                                       `9` = "CD8pos_T_cell",
                                       `10` = "CD4pos_Tconv",
                                       `11` = "CD4pos_Tconv",
                                       `12` = "CD4pos_Treg",
                                       `13` = "CD8pos_T_cell"
                                                                         )
               )

    meta$refined_lymph_cell_annot <- as.character(meta$refined_lymph_cell_annot)

    table(meta$refined_lymph_cell_annot,useNA="always")

    ## 
    ##  CD4pos_Tconv   CD4pos_Treg CD8pos_T_cell       NK_cell          <NA> 
    ##          1708            49          1417           183             0

    ser_comb@meta.data$refined_nk_t_cell_annot <- meta$refined_lymph_cell_annot

    DimPlot(ser_comb,group.by="refined_nk_t_cell_annot")

![](04_refined_cell_annotations-exported_files/figure-markdown_strict/unnamed-chunk-6-1.png)

    refined_nk_t_results <- ser_comb

Re-load datasets
----------------

    ser_comb <- readRDS("../../03_output/SARC028_CosMx_and_scRNAseq_object.rds")

    DimPlot(ser_comb,group.by="initial_cell_annot",label=T)

![](04_refined_cell_annotations-exported_files/figure-markdown_strict/unnamed-chunk-7-1.png)

Re-cluster Plasmablasts and pDCs
--------------------------------

We note here that we do not need to re-integrate plasmablasts and pDCs
to separate them by cell types.

    ser_comb <- ser_comb[,grepl("Plasmablasts/pDCs",ser_comb@meta.data$initial_cell_annot)]

    table(ser_comb@meta.data$initial_cell_annot)

    ## 
    ##       Tumor_cells_c0       Tumor_cells_c1       Tumor_cells_c2 
    ##                    0                    0                    0 
    ##       Tumor_cells_c3          Osteoclasts       Tumor_cells_c4 
    ##                    0                    0                    0 
    ##       T_and_NK_cells       Tumor_cells_c5       Tumor_cells_c6 
    ##                    0                    0                    0 
    ## CD14posCD16pos_macro       Tumor_cells_c7            Monocytes 
    ##                    0                    0                    0 
    ##             CD1C_DCs    Plasmablasts/pDCs          Fibroblasts 
    ##                    0                 1045                    0 
    ##              B_cells   C1QCpos_macrophage 
    ##                    0                    0

    DimPlot(ser_comb,group.by="initial_cell_annot")

![](04_refined_cell_annotations-exported_files/figure-markdown_strict/unnamed-chunk-8-1.png)

    DimPlot(ser_comb,group.by="data_type")

![](04_refined_cell_annotations-exported_files/figure-markdown_strict/unnamed-chunk-8-2.png)

    FeaturePlot(ser_comb,c("MZB1","IL3RA","JCHAIN"))

![](04_refined_cell_annotations-exported_files/figure-markdown_strict/unnamed-chunk-8-3.png)

    DefaultAssay(ser_comb) <- "RNA"

    ser_comb <- ser_comb %>%
        NormalizeData(.) %>%
        FindVariableFeatures(.) %>%
        ScaleData(.) %>%
        RunPCA(., npcs = 30, verbose = FALSE)

    ## Centering and scaling data matrix

    ElbowPlot(ser_comb)

![](04_refined_cell_annotations-exported_files/figure-markdown_strict/unnamed-chunk-8-4.png)

    ser_comb <- RunUMAP(ser_comb, reduction = "pca", dims = 1:15)

    ## 10:51:59 UMAP embedding parameters a = 0.9922 b = 1.112

    ## 10:51:59 Read 1045 rows and found 15 numeric columns

    ## 10:51:59 Using Annoy for neighbor search, n_neighbors = 30

    ## 10:51:59 Building Annoy index with metric = cosine, n_trees = 50

    ## 0%   10   20   30   40   50   60   70   80   90   100%

    ## [----|----|----|----|----|----|----|----|----|----|

    ## **************************************************|
    ## 10:51:59 Writing NN index file to temp file /tmp/RtmpjvQNBv/filee6ed1d00ba5c
    ## 10:51:59 Searching Annoy index using 1 thread, search_k = 3000
    ## 10:51:59 Annoy recall = 100%
    ## 10:51:59 Commencing smooth kNN distance calibration using 1 thread with target n_neighbors = 30
    ## 10:52:00 Initializing from normalized Laplacian + noise (using irlba)
    ## 10:52:00 Commencing optimization for 500 epochs, with 39866 positive edges
    ## 10:52:05 Optimization finished

    ser_comb <- FindNeighbors(ser_comb, reduction = "pca", dims = 1:15)

    ## Computing nearest neighbor graph
    ## Computing SNN

    ser_comb <- FindClusters(ser_comb, resolution = c(0.3,0.5,0.7,1))

    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 1045
    ## Number of edges: 30845
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.9334
    ## Number of communities: 9
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 1045
    ## Number of edges: 30845
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.9114
    ## Number of communities: 11
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 1045
    ## Number of edges: 30845
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8905
    ## Number of communities: 12
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 1045
    ## Number of edges: 30845
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8615
    ## Number of communities: 13
    ## Elapsed time: 0 seconds

Refine Plasmablasts and pDCs
----------------------------

We note that clusters 0 and 5 do not robuslty express markers of either
plasmablasts or pDCs and are therefore annotated as NA.

    DimPlot(ser_comb,group.by="RNA_snn_res.0.3",label=T)

![](04_refined_cell_annotations-exported_files/figure-markdown_strict/unnamed-chunk-9-1.png)

    Idents(ser_comb) <- "RNA_snn_res.0.3"

    degs <- FindAllMarkers(ser_comb)

    ## Calculating cluster 0

    ## Calculating cluster 1

    ## Calculating cluster 2

    ## Calculating cluster 3

    ## Calculating cluster 4

    ## Calculating cluster 5

    ## Calculating cluster 6

    ## Calculating cluster 7

    ## Calculating cluster 8

    degs %>%
        filter(p_val_adj<0.05) %>%
        filter(avg_log2FC>0) %>%
        split(.$cluster) %>%
        lapply(.,function(x) x %>% arrange(desc(avg_log2FC)) %>% head(n=10))

    ## $`0`
    ##                 p_val avg_log2FC pct.1 pct.2     p_val_adj cluster    gene
    ## COL18A1 9.427756e-115   3.258422 0.781 0.059 9.003507e-112       0 COL18A1
    ## CDH1    1.500045e-100   3.110770 0.651 0.036  1.432543e-97       0    CDH1
    ## SLC40A1 8.830307e-124   2.979613 0.912 0.100 8.432944e-121       0 SLC40A1
    ## CPB1     1.051307e-54   2.899156 0.302 0.005  1.003998e-51       0    CPB1
    ## EPOR    3.444688e-107   2.772499 0.842 0.093 3.289677e-104       0    EPOR
    ## IL7      3.390707e-66   2.689915 0.409 0.016  3.238125e-63       0     IL7
    ## CLU      8.674630e-43   2.646255 0.284 0.014  8.284272e-40       0     CLU
    ## CSF1     2.983725e-63   2.645829 0.428 0.024  2.849458e-60       0    CSF1
    ## CYTOR   6.079872e-113   2.607196 0.953 0.141 5.806278e-110       0   CYTOR
    ## CD36     1.180835e-90   2.597652 0.898 0.175  1.127697e-87       0    CD36
    ## 
    ## $`1`
    ##                p_val avg_log2FC pct.1 pct.2    p_val_adj cluster  gene
    ## BMP3    8.293093e-57   4.500163 0.326 0.009 7.919904e-54       1  BMP3
    ## ROR1.1  6.751701e-60   3.613809 0.605 0.091 6.447874e-57       1  ROR1
    ## HCK.1   4.044812e-63   3.522827 0.605 0.084 3.862795e-60       1   HCK
    ## FCRLA.1 9.488956e-41   3.043521 0.519 0.105 9.061953e-38       1 FCRLA
    ## TCL1A.1 3.599459e-77   2.757877 0.853 0.163 3.437483e-74       1 TCL1A
    ## CD19.1  4.506106e-44   2.456709 0.674 0.181 4.303331e-41       1  CD19
    ## P2RX5.1 4.598355e-19   2.350052 0.481 0.198 4.391429e-16       1 P2RX5
    ## IRF4.1  3.194050e-39   2.057706 0.752 0.267 3.050318e-36       1  IRF4
    ## MAML2.1 1.387612e-17   1.985493 0.496 0.209 1.325169e-14       1 MAML2
    ## IL4R    1.951958e-05   1.935634 0.155 0.059 1.864120e-02       1  IL4R
    ## 
    ## $`2`
    ##               p_val avg_log2FC pct.1 pct.2    p_val_adj cluster   gene
    ## MPO    3.743518e-44   4.489211 0.370 0.031 3.575059e-41       2    MPO
    ## ELANE  1.274915e-26   4.273840 0.210 0.015 1.217544e-23       2  ELANE
    ## CSF3R  2.746060e-52   3.579373 0.387 0.024 2.622487e-49       2  CSF3R
    ## CTSG   1.602242e-26   3.517140 0.202 0.013 1.530141e-23       2   CTSG
    ## AZU1   5.596005e-13   2.977838 0.244 0.060 5.344185e-10       2   AZU1
    ## CST7   6.299776e-12   2.743065 0.210 0.049 6.016286e-09       2   CST7
    ## BCL2.1 4.494531e-41   2.588736 0.639 0.148 4.292277e-38       2   BCL2
    ## SUCNR1 1.759224e-24   2.431785 0.151 0.005 1.680059e-21       2 SUCNR1
    ## CFD    3.164676e-08   2.409144 0.210 0.070 3.022266e-05       2    CFD
    ## HGF    9.363169e-19   2.358435 0.109 0.003 8.941827e-16       2    HGF
    ## 
    ## $`3`
    ##                  p_val avg_log2FC pct.1 pct.2     p_val_adj cluster    gene
    ## DNTT.3   6.580062e-118   3.618781 0.966 0.115 6.283959e-115       3    DNTT
    ## COL5A1.3  4.119209e-36   2.951988 0.441 0.071  3.933844e-33       3  COL5A1
    ## CD34.2    2.881521e-24   2.724378 0.347 0.067  2.751852e-21       3    CD34
    ## CD9.3     1.011896e-67   2.579699 0.898 0.248  9.663603e-65       3     CD9
    ## CLEC14A   5.827253e-14   2.562113 0.102 0.006  5.565026e-11       3 CLEC14A
    ## MT1X.1    2.092011e-20   2.470209 0.500 0.196  1.997870e-17       3    MT1X
    ## VPREB3.2  6.474043e-74   2.224056 0.983 0.289  6.182711e-71       3  VPREB3
    ## ADGRE2.2  1.104237e-15   2.179680 0.347 0.104  1.054546e-12       3  ADGRE2
    ## LTB.3     4.002162e-39   2.058552 0.746 0.216  3.822064e-36       3     LTB
    ## PECAM1.3  1.425915e-29   2.057709 0.619 0.211  1.361749e-26       3  PECAM1
    ## 
    ## $`4`
    ##                   p_val avg_log2FC pct.1 pct.2    p_val_adj cluster       gene
    ## TCL1A.4    4.322755e-69   2.387362 0.956 0.161 4.128231e-66       4      TCL1A
    ## IL7R.2     1.543920e-43   2.248612 0.482 0.056 1.474444e-40       4       IL7R
    ## MS4A1.4    1.104678e-37   2.241337 0.421 0.048 1.054967e-34       4      MS4A1
    ## CRIP1.1    2.687176e-62   2.144053 0.974 0.301 2.566253e-59       4      CRIP1
    ## ROR1.4     5.850578e-60   1.880571 0.719 0.085 5.587302e-57       4       ROR1
    ## FCRLA.4    4.041761e-38   1.789132 0.596 0.102 3.859882e-35       4      FCRLA
    ## CD19.2     8.552538e-45   1.763248 0.825 0.171 8.167674e-42       4       CD19
    ## ST6GALNAC3 1.393307e-20   1.653332 0.325 0.059 1.330608e-17       4 ST6GALNAC3
    ## CD24.3     2.607612e-42   1.652494 0.974 0.270 2.490270e-39       4       CD24
    ## CD79A.3    3.997146e-36   1.497128 1.000 0.344 3.817274e-33       4      CD79A
    ## 
    ## $`5`
    ##                 p_val avg_log2FC pct.1 pct.2    p_val_adj cluster   gene
    ## HBA1.5   2.400411e-36   2.484797 0.728 0.240 2.292393e-33       5   HBA1
    ## LMNA.5   2.263913e-47   2.409509 0.835 0.289 2.162037e-44       5   LMNA
    ## EPCAM.5  1.469654e-16   2.380877 0.427 0.143 1.403520e-13       5  EPCAM
    ## LGALS3.5 2.271892e-18   2.329432 0.476 0.163 2.169657e-15       5 LGALS3
    ## IFI27.5  8.444772e-10   2.312057 0.272 0.091 8.064757e-07       5  IFI27
    ## TPM1.5   4.037292e-24   2.271007 0.583 0.220 3.855614e-21       5   TPM1
    ## KRT13.5  3.588591e-06   2.102830 0.214 0.085 3.427104e-03       5  KRT13
    ## S100A6.5 1.663677e-60   2.028485 0.990 0.534 1.588811e-57       5 S100A6
    ## S100A4.5 2.280060e-57   1.948190 0.990 0.538 2.177457e-54       5 S100A4
    ## MT2A.4   5.813149e-14   1.931911 0.505 0.243 5.551557e-11       5   MT2A
    ## 
    ## $`6`
    ##                  p_val avg_log2FC pct.1 pct.2     p_val_adj cluster   gene
    ## CD68.5    3.175284e-97   3.708823 0.788 0.061  3.032396e-94       6   CD68
    ## GZMB.6   9.474205e-110   3.505345 0.859 0.052 9.047866e-107       6   GZMB
    ## IL18R1    5.618511e-77   3.439523 0.465 0.013  5.365678e-74       6 IL18R1
    ## FCER1G.5 1.664268e-102   3.325730 0.909 0.081  1.589376e-99       6 FCER1G
    ## CMKLR1    4.122786e-62   3.265274 0.364 0.008  3.937261e-59       6 CMKLR1
    ## IL3RA.5   1.180379e-97   3.222962 0.929 0.101  1.127262e-94       6  IL3RA
    ## TYROBP.5  1.293931e-91   3.134055 0.909 0.097  1.235704e-88       6 TYROBP
    ## TPM2.4    6.300545e-81   3.108192 0.747 0.067  6.017021e-78       6   TPM2
    ## DPP4      2.942054e-42   3.029211 0.242 0.005  2.809662e-39       6   DPP4
    ## CSF2RA.5  3.357461e-84   2.873070 0.707 0.049  3.206375e-81       6 CSF2RA
    ## 
    ## $`7`
    ##                  p_val avg_log2FC pct.1 pct.2    p_val_adj cluster    gene
    ## GADD45B.6 4.934468e-09  1.2337798 0.679 0.644 4.712417e-06       7 GADD45B
    ## ZFP36.5   4.989911e-07  1.0378262 0.628 0.631 4.765365e-04       7   ZFP36
    ## MALAT1.3  7.551065e-30  0.8994265 0.949 0.923 7.211267e-27       7  MALAT1
    ## FOS.2     1.238438e-10  0.5731047 0.808 0.871 1.182708e-07       7     FOS
    ## 
    ## $`8`
    ##                   p_val avg_log2FC pct.1 pct.2    p_val_adj cluster      gene
    ## IGHA1      6.130390e-65   5.915069 0.443 0.015 5.854523e-62       8     IGHA1
    ## IGHG2      3.040081e-57   5.604191 0.271 0.001 2.903278e-54       8     IGHG2
    ## IGHG1      6.737664e-49   5.423609 0.300 0.007 6.434469e-46       8     IGHG1
    ## CD27.1     1.037738e-38   3.601206 0.429 0.039 9.910403e-36       8      CD27
    ## GPR183.6   2.011753e-37   3.572967 0.543 0.076 1.921224e-34       8    GPR183
    ## WNT5B      1.757243e-40   3.567446 0.314 0.014 1.678167e-37       8     WNT5B
    ## JCHAIN.5   2.193488e-56   3.535042 0.900 0.208 2.094781e-53       8    JCHAIN
    ## TNFRSF13B  2.685489e-43   3.495591 0.314 0.012 2.564642e-40       8 TNFRSF13B
    ## TNFRSF17.2 7.180189e-42   3.296621 0.500 0.051 6.857080e-39       8  TNFRSF17
    ## PNOC.1     1.749214e-34   3.270943 0.329 0.023 1.670500e-31       8      PNOC

    FeaturePlot(ser_comb,c("MS4A1","IL3RA"))

![](04_refined_cell_annotations-exported_files/figure-markdown_strict/unnamed-chunk-9-2.png)

    b_pdc_meta <- ser_comb@meta.data %>%
        mutate(refined_b_pdc_annot=ifelse(RNA_snn_res.0.3 %in% c(4,1,7,3),"Plasmablasts","pDCs")) %>%
        mutate(refined_b_pdc_annot=ifelse(RNA_snn_res.0.3 %in% c(0,5),NA,refined_b_pdc_annot))

    table(b_pdc_meta$refined_b_pdc_annot,useNA="always")

    ## 
    ##         pDCs Plasmablasts         <NA> 
    ##          288          439          318

Add refined annotations to original dataset
-------------------------------------------

    ## Reload original data
    ser_comb <- readRDS("../../03_output/SARC028_CosMx_and_scRNAseq_object.rds")

    meta_orig <- ser_comb@meta.data %>%
        as_tibble(.,rownames="cb")

    meta_new1 <- refined_nk_t_results@meta.data %>%
        as_tibble(.,rownames="cb") %>%
        select(cb,refined_nk_t_cell_annot)

    meta_new2 <- b_pdc_meta %>%
        as_tibble(.,rownames="cb") %>%
        select(cb,refined_b_pdc_annot)

    meta_all <- left_join(meta_orig,meta_new1,by="cb") %>%
        left_join(.,meta_new2,by="cb") %>%
        mutate(final_cell_annot=ifelse(is.na(refined_nk_t_cell_annot),
                                               as.character(initial_cell_annot),
                                               refined_nk_t_cell_annot)
               ) %>%
        mutate(final_cell_annot=ifelse(is.na(refined_b_pdc_annot),
                                               as.character(final_cell_annot),
                                               refined_b_pdc_annot)
               )

    table(meta_all$final_cell_annot,useNA="always")

    ## 
    ##              B_cells   C1QCpos_macrophage CD14posCD16pos_macro 
    ##                  324                  258                 2320 
    ##             CD1C_DCs         CD4pos_Tconv          CD4pos_Treg 
    ##                 1230                 1708                   49 
    ##        CD8pos_T_cell          Fibroblasts            Monocytes 
    ##                 1417                  435                 1417 
    ##              NK_cell          Osteoclasts                 pDCs 
    ##                  183                 3846                  288 
    ##         Plasmablasts    Plasmablasts/pDCs       Tumor_cells_c0 
    ##                  439                  318                12992 
    ##       Tumor_cells_c1       Tumor_cells_c2       Tumor_cells_c3 
    ##                12940                 7758                 4549 
    ##       Tumor_cells_c4       Tumor_cells_c5       Tumor_cells_c6 
    ##                 3512                 2798                 2479 
    ##       Tumor_cells_c7                 <NA> 
    ##                 2268                    0

    ser_comb@meta.data$final_cell_annot <- meta_all$final_cell_annot

    ser_comb <- ser_comb[,!ser_comb@meta.data$final_cell_annot=="Plasmablasts/pDCs"]

    DimPlot(ser_comb,group.by="final_cell_annot")

![](04_refined_cell_annotations-exported_files/figure-markdown_strict/unnamed-chunk-10-1.png)

Save results and session info
-----------------------------

    # saveRDS(ser_comb,file="../../03_output/SARC028_CosMx_and_scRNAseq_object_refined_annotations.rds")

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
    ##  [61] promises_1.2.0.1       spatstat.utils_3.0-1   parallel_4.2.0        
    ##  [64] RColorBrewer_1.1-3     yaml_2.3.6             reticulate_1.26       
    ##  [67] pbapply_1.6-0          gridExtra_2.3          stringi_1.7.8         
    ##  [70] highr_0.9              rlang_1.1.1            pkgconfig_2.0.3       
    ##  [73] matrixStats_0.63.0     evaluate_0.15          lattice_0.20-45       
    ##  [76] ROCR_1.0-11            tensor_1.5             labeling_0.4.2        
    ##  [79] htmlwidgets_1.5.4      cowplot_1.1.1          tidyselect_1.2.0      
    ##  [82] parallelly_1.32.1      RcppAnnoy_0.0.20       plyr_1.8.8            
    ##  [85] magrittr_2.0.3         R6_2.5.1               generics_0.1.3        
    ##  [88] pillar_1.8.1           withr_2.5.0            fitdistrplus_1.1-8    
    ##  [91] survival_3.5-3         abind_1.4-5            sp_1.5-1              
    ##  [94] future.apply_1.10.0    KernSmooth_2.23-20     utf8_1.2.2            
    ##  [97] spatstat.geom_3.0-3    plotly_4.10.1          tzdb_0.3.0            
    ## [100] grid_4.2.0             data.table_1.14.6      digest_0.6.30         
    ## [103] xtable_1.8-4           httpuv_1.6.6           munsell_0.5.0         
    ## [106] viridisLite_0.4.1
