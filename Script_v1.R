##Packages----

#Installation
install.packages("Seurat")

#Load Packages
library(Seurat)




#Load Seurat Object----
seu <- readRDS("seurat_ass4.rds")

seu
class(seu)


#Inspect Structure

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

#There are missing mouse_id values spread throughout the data set. Therefore, later on, we will use biosample_id for DE. For now, we may continue with the regular intended pipeline.


seu <- FindVariableFeatures(seu)
length(VariableFeatures(seu))

seu <- ScaleData(seu, features = VariableFeatures(seu), verbose = FALSE)
saveRDS(seu, file = "seu_after_scaling.rds", compress = TRUE)

seu <- RunPCA(seu, features = VariableFeatures(seu), verbose = FALSE)

ElbowPlot(seu, ndims = 30)
#1:12 seems to be the most fitting

seu <- FindNeighbors(seu, dims = 1:12, verbose = FALSE)

seu <- FindClusters(seu, resolution = 0.3, verbose = FALSE)

seu <- RunUMAP(seu, dims = 1:12, verbose = FALSE)


DimPlot(seu, reduction = "umap", label = TRUE)
table(Idents(seu))


DimPlot(seu, reduction = "umap", group.by = "organ_custom")

DimPlot(seu, reduction = "umap", group.by = "time")
