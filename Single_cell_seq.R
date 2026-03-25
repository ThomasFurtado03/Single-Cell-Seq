##Packages----

#Install (if needed)
install.packages("Seurat")

#Load Packages
library(Seurat)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)

##Load Seurat Object----

seu <- readRDS("seurat_ass4.rds")

#Load pre-scaled object (to prevent having to re-scale every session) (Overwrites previous)
seu <- readRDS("seu_after_scaling.rds")


seu
class(seu)


##Data Inspection and QC Checks----

Assays(seu)

colnames(seu@meta.data)
head(seu@meta.data)

table(Idents(seu))

dim(seu[["RNA"]]@counts)
dim(seu[["RNA"]]@data)

table(seu$organ_custom)
table(seu$time)
table(seu$biosample_id)
table(seu$orig.ident)

#Check for NA
table(seu$mouse_id, useNA = "ifany")
sum(seu$mouse_id == "" | is.na(seu$mouse_id))

table(seu$mouse_id == "", seu$time)
table(seu$mouse_id == "", seu$organ_custom)
table(seu$mouse_id == "", seu$biosample_id)
table(seu$mouse_id == "", seu$orig.ident)

head(seu@meta.data[seu$mouse_id == "", ])

#Note many missing mouse_id, may need to use biosample_id for DE


##Feature Selection and Scaling----

seu <- FindVariableFeatures(seu)
length(VariableFeatures(seu))

seu <- ScaleData(seu, features = VariableFeatures(seu), verbose = FALSE)

#Save scaled object
saveRDS(seu, file = "seu_after_scaling.rds", compress = TRUE)



##Dimensionality Reduction and Clustering----

seu <- RunPCA(seu, features = VariableFeatures(seu), verbose = FALSE)

ElbowPlot(seu, ndims = 30)
#1:12 seems to be the most fitting

seu <- FindNeighbors(seu, dims = 1:12, verbose = FALSE)
seu <- FindClusters(seu, resolution = 0.3, verbose = FALSE)

seu <- RunUMAP(seu, dims = 1:12, verbose = FALSE)



##UMAP Visualization and Cluster Exploration----

DimPlot(seu, reduction = "umap", label = TRUE)
table(Idents(seu))


DimPlot(seu, reduction = "umap", group.by = "organ_custom")
DimPlot(seu, reduction = "umap", group.by = "time")



#Lets choose a reasonably sized cluster that also shows some time variation.
#To do this lets look at how each cluster is distributed over time.
table(Idents(seu), seu$time)

#Cluster 5 has quite a lot of variation across the studied time periods. I will choose to investigate this cluster.



##Differential Expression: Cluster 5 -> D05 vs D14----

cluster5 <- subset(seu, idents = 5)
table(cluster5$time)


Idents(cluster5) <- cluster5$time

markers_5 <- FindMarkers(
  cluster5,
  ident.1 = "D14",
  ident.2 = "D05",
  min.pct = 0.1,
  logfc.threshold = 0.25
)

head(markers_5)

#Obp1a has very low (-310) avg_log2FC, meaning it is much higher in D05 than D14

#Most of the top hits are "Rps___ / Rpl___" ribosomal genes, which are non-specific and likely dominate other meaningful results. Lets filter for biologically meaningful genes.

markers_5_filtered <- markers_5[
  markers_5$p_val_adj < 0.05 &
    abs(markers_5$avg_log2FC) > 0.5,
]


markers_5_clean <- markers_5_filtered[
  !grepl("^Rps|^Rpl", rownames(markers_5_filtered)),
]

head(markers_5_clean)
write.csv(markers_5_clean, "DE_results.csv")

#Almost all of the top genes have negative log2FC values, meaning that they are higher in D05, than in D14.

#Instinctively, cluster5 genes at D05 may be expressing immune-related genes, which come to be reduced come D14.


##Visualize select genes----

select_genes <- c("Ifitm2", "Ifi27l2a", "Il7r", "Obp1a")
FeaturePlot(cluster5, features = select_genes, min.cutoff = "q10", max.cutoff = "q90")


##Functional Enrichment (ORA)----

sig_genes <- rownames(markers_5_clean)

ego <- enrichGO(
  gene = sig_genes,
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "bp"
)


#Visualize
dotplot(ego, showCategory = 9) +
  theme(
    plot.margin = margin(t = 0, r = 5, b = 0, l = 0),
    axis.text.y = element_text(size = 7, lineheight = 0.7)
  )

##Annotation of Cluster 5----

#Due to computational strain, annotation will be done using a small subset of the data while still commiting a true marker-based annotation step.
set.seed(123)

Idents(seu) <- "seurat_clusters"

cells_use <- c(
  sample(WhichCells(seu, idents = 5), 1500),
  sample(WhichCells(seu, idents = setdiff(levels(Idents(seu)), "5")), 6000)
)

seu_small <- subset(seu, cells = cells_use)

cluster5_markers <- FindMarkers(
  seu_small,
  ident.1 = "5",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

cluster5_markers_clean <- cluster5_markers[
  !grepl("^Rps|^Rpl|^mt-", rownames(cluster5_markers)),
]

head(cluster5_markers_clean, 15)

top15 <- cluster5_markers_clean[order(-cluster5_markers_clean$avg_log2FC), ][1:15, ]
top15_df <- data.frame(
  Gene = rownames(top15),
  log2FC = round(top15$avg_log2FC, 2),
  pct_D14 = round(top15$pct.1, 2),
  pct_D05 = round(top15$pct.2, 2)
)

top15_df
write.csv(top15_df, "cluster5_top15_markers.csv")

#Genes such as Cd79a, Cd79b, Ms4a1, Pax5, Igkc, Ighm, etc. are all classic B-cell markers!
#Therefore Cluster 5 is a B-cell population



##Biological Interpretation and Summary----

#To summarize: 
#From Markers -> B-cell genes identified, Cluster 5 -> B-cells

#From DE -> Interferon/anti-viral genes identified -> active response early (D05)

#From ORA -> Immune activation, antigen presentation, viral response


#Likely conclusion:
#Cluster 5 very likely represents a B-cell population which exhibits strong antiviral immune response, which is high during early infection (D05) and depletes by later stages (D14)


