library(dplyr)
library(Seurat)
library(future)
plan("multisession", workers =5)
#plan("sequential")
library(ggplot2)
options(future.seed = TRUE)
options(future.globals.maxSize = 100 * 1024^3)

path.1 <- "/rsrch3/home/genomic_med/ssatpati/LUNG_Mixed_Histo/20251016__203455__LUNG_X5K_9MJMHF_Batch1_10162025/output-XETG00458__0067534__1__20251016__203659/"
#xenium.obj.1 <- LoadXenium(path.1, fov = "fov")
xenium.obj.1 <- ReadXenium(data.dir = path.1,  outs = c("segmentation_method", "matrix", "microns"),  type = "centroids",  mols.qv.threshold = 20)
gene_expr1 <- xenium.obj.1$matrix[["Gene Expression"]]
xenium.obj.1 <- CreateSeuratObject(counts = gene_expr1, assay = "Xenium")
xenium.obj.1 <- subset(xenium.obj.1, subset = nCount_Xenium > 0)

#path.2 <- "/rsrch3/home/genomic_med/ssatpati/LUNG_Mixed_Histo/20251016__203455__LUNG_X5K_9MJMHF_Batch1_10162025/output-XETG00458__0067534__2__20251016__203659/"
#xenium.obj.2 <- LoadXenium(path.2, fov = "fov")
#xenium.obj.2 <- subset(xenium.obj.2, subset = nCount_Xenium > 0)

path.3 <- "/rsrch3/home/genomic_med/ssatpati/LUNG_Mixed_Histo/20251016__203455__LUNG_X5K_9MJMHF_Batch1_10162025/output-XETG00458__0067534__3__20251016__203659/"
#xenium.obj.3 <- LoadXenium(path.3, fov = "fov")
xenium.obj.3 <- ReadXenium(data.dir = path.3,  outs = c("segmentation_method", "matrix", "microns"),  type = "centroids",  mols.qv.threshold = 20)
gene_expr3 <- xenium.obj.3$matrix[["Gene Expression"]]
xenium.obj.3 <- CreateSeuratObject(counts = gene_expr3, assay = "Xenium")
xenium.obj.3 <- subset(xenium.obj.3, subset = nCount_Xenium > 0)

path.4 <- "/rsrch3/home/genomic_med/ssatpati/LUNG_Mixed_Histo/20251016__203455__LUNG_X5K_9MJMHF_Batch1_10162025/output-XETG00458__0067534__4__20251016__203659/"
#xenium.obj.4 <- LoadXenium(path.4, fov = "fov")
xenium.obj.4 <- ReadXenium(data.dir = path.4,  outs = c("segmentation_method", "matrix", "microns"),  type = "centroids",  mols.qv.threshold = 20)
gene_expr4 <- xenium.obj.4$matrix[["Gene Expression"]]
xenium.obj.4 <- CreateSeuratObject(counts = gene_expr4, assay = "Xenium")
xenium.obj.4 <- subset(xenium.obj.4, subset = nCount_Xenium > 0)

path.5 <- "/rsrch3/home/genomic_med/ssatpati/LUNG_Mixed_Histo/20251016__203455__LUNG_X5K_9MJMHF_Batch1_10162025/output-XETG00458__0067555__1__20251016__203659/"
#xenium.obj.5 <- LoadXenium(path.5, fov = "fov")
xenium.obj.5 <- ReadXenium(data.dir = path.5,  outs = c("segmentation_method", "matrix", "microns"),  type = "centroids",  mols.qv.threshold = 20)
gene_expr5 <- xenium.obj.5$matrix[["Gene Expression"]]
xenium.obj.5 <- CreateSeuratObject(counts = gene_expr5, assay = "Xenium")
xenium.obj.5 <- subset(xenium.obj.5, subset = nCount_Xenium > 0)

path.6 <- "/rsrch3/home/genomic_med/ssatpati/LUNG_Mixed_Histo/20251016__203455__LUNG_X5K_9MJMHF_Batch1_10162025/output-XETG00458__0067555__2__20251016__203659/"
#xenium.obj.6 <- LoadXenium(path.6, fov = "fov")
xenium.obj.6 <- ReadXenium(data.dir = path.6,  outs = c("segmentation_method", "matrix", "microns"),  type = "centroids",  mols.qv.threshold = 20)
gene_expr6 <- xenium.obj.6$matrix[["Gene Expression"]]
xenium.obj.6 <- CreateSeuratObject(counts = gene_expr6, assay = "Xenium")
xenium.obj.6 <- subset(xenium.obj.6, subset = nCount_Xenium > 0)

path.7 <- "/rsrch3/home/genomic_med/ssatpati/LUNG_Mixed_Histo/20251016__203455__LUNG_X5K_9MJMHF_Batch1_10162025/output-XETG00458__0067555__3__20251016__203659/"
#xenium.obj.7 <- LoadXenium(path.7, fov = "fov")
xenium.obj.7 <- ReadXenium(data.dir = path.7,  outs = c("segmentation_method", "matrix", "microns"),  type = "centroids",  mols.qv.threshold = 20)
gene_expr7 <- xenium.obj.7$matrix[["Gene Expression"]]
xenium.obj.7 <- CreateSeuratObject(counts = gene_expr7, assay = "Xenium")
xenium.obj.7 <- subset(xenium.obj.7, subset = nCount_Xenium > 0)

path.8 <- "/rsrch3/home/genomic_med/ssatpati/LUNG_Mixed_Histo/20251016__203455__LUNG_X5K_9MJMHF_Batch1_10162025/output-XETG00458__0067555__4__20251016__203659/"
#xenium.obj.8 <- LoadXenium(path.8, fov = "fov")
xenium.obj.8 <- ReadXenium(data.dir = path.8,  outs = c("segmentation_method", "matrix", "microns"),  type = "centroids",  mols.qv.threshold = 20)
gene_expr8 <- xenium.obj.8$matrix[["Gene Expression"]]
xenium.obj.8 <- CreateSeuratObject(counts = gene_expr8, assay = "Xenium")
xenium.obj.8 <- subset(xenium.obj.8, subset = nCount_Xenium > 0)

path.9 <- "/rsrch3/home/genomic_med/ssatpati/LUNG_Mixed_Histo/20251029__211811__X5K_9MJMHF_Batch2_Lung_10292025/output-XETG00458__0073238__1__20251029__212017/"
#xenium.obj.9 <- LoadXenium(path.9, fov = "fov")
xenium.obj.9 <- ReadXenium(data.dir = path.9,  outs = c("segmentation_method", "matrix", "microns"),  type = "centroids",  mols.qv.threshold = 20)
gene_expr9 <- xenium.obj.9$matrix[["Gene Expression"]]
xenium.obj.9 <- CreateSeuratObject(counts = gene_expr9, assay = "Xenium")
xenium.obj.9 <- subset(xenium.obj.9, subset = nCount_Xenium > 0)

path.10 <- "/rsrch3/home/genomic_med/ssatpati/LUNG_Mixed_Histo/20251029__211811__X5K_9MJMHF_Batch2_Lung_10292025/output-XETG00458__0073238__2__20251029__212017/"
#xenium.obj.10 <- LoadXenium(path.10, fov = "fov")
xenium.obj.10 <- ReadXenium(data.dir = path.10,  outs = c("segmentation_method", "matrix", "microns"),  type = "centroids",  mols.qv.threshold = 20)
gene_expr10 <- xenium.obj.10$matrix[["Gene Expression"]]
xenium.obj.10 <- CreateSeuratObject(counts = gene_expr10, assay = "Xenium")
xenium.obj.10 <- subset(xenium.obj.10, subset = nCount_Xenium > 0)

path.11 <- "/rsrch3/home/genomic_med/ssatpati/LUNG_Mixed_Histo/20251029__211811__X5K_9MJMHF_Batch2_Lung_10292025/output-XETG00458__0073238__3__20251029__212017/"
#xenium.obj.11 <- LoadXenium(path.11, fov = "fov")
xenium.obj.11 <- ReadXenium(data.dir = path.11,  outs = c("segmentation_method", "matrix", "microns"),  type = "centroids",  mols.qv.threshold = 20)
gene_expr11 <- xenium.obj.11$matrix[["Gene Expression"]]
xenium.obj.11 <- CreateSeuratObject(counts = gene_expr11, assay = "Xenium")
xenium.obj.11 <- subset(xenium.obj.11, subset = nCount_Xenium > 0)

path.12 <- "/rsrch3/home/genomic_med/ssatpati/LUNG_Mixed_Histo/20251029__211811__X5K_9MJMHF_Batch2_Lung_10292025/output-XETG00458__0073238__4__20251029__212017/"
#xenium.obj.12 <- LoadXenium(path.12, fov = "fov")
xenium.obj.12 <- ReadXenium(data.dir = path.12,  outs = c("segmentation_method", "matrix", "microns"),  type = "centroids",  mols.qv.threshold = 20)
gene_expr12 <- xenium.obj.12$matrix[["Gene Expression"]]
xenium.obj.12 <- CreateSeuratObject(counts = gene_expr12, assay = "Xenium")
xenium.obj.12 <- subset(xenium.obj.12, subset = nCount_Xenium > 0)

path.13 <- "/rsrch3/home/genomic_med/ssatpati/LUNG_Mixed_Histo/20251029__211811__X5K_9MJMHF_Batch2_Lung_10292025/output-XETG00458__0073241__1__20251029__212017/"
#xenium.obj.13 <- LoadXenium(path.13, fov = "fov")
xenium.obj.13 <- ReadXenium(data.dir = path.13,  outs = c("segmentation_method", "matrix", "microns"),  type = "centroids",  mols.qv.threshold = 20)
gene_expr13 <- xenium.obj.13$matrix[["Gene Expression"]]
xenium.obj.13 <- CreateSeuratObject(counts = gene_expr13, assay = "Xenium")
xenium.obj.13 <- subset(xenium.obj.13, subset = nCount_Xenium > 0)

path.14 <- "/rsrch3/home/genomic_med/ssatpati/LUNG_Mixed_Histo/20251029__211811__X5K_9MJMHF_Batch2_Lung_10292025/output-XETG00458__0073241__2__20251029__212017/"
#xenium.obj.14 <- LoadXenium(path.14, fov = "fov")
xenium.obj.14 <- ReadXenium(data.dir = path.14,  outs = c("segmentation_method", "matrix", "microns"),  type = "centroids",  mols.qv.threshold = 20)
gene_expr14 <- xenium.obj.14$matrix[["Gene Expression"]]
xenium.obj.14 <- CreateSeuratObject(counts = gene_expr14, assay = "Xenium")
xenium.obj.14 <- subset(xenium.obj.14, subset = nCount_Xenium > 0)

path.15 <- "/rsrch3/home/genomic_med/ssatpati/LUNG_Mixed_Histo/20251029__211811__X5K_9MJMHF_Batch2_Lung_10292025/output-XETG00458__0073241__3__20251029__212017/"
#xenium.obj.15 <- LoadXenium(path.15, fov = "fov")
xenium.obj.15 <- ReadXenium(data.dir = path.15,  outs = c("segmentation_method", "matrix", "microns"),  type = "centroids",  mols.qv.threshold = 20)
gene_expr15 <- xenium.obj.15$matrix[["Gene Expression"]]
xenium.obj.15 <- CreateSeuratObject(counts = gene_expr15, assay = "Xenium")
xenium.obj.15 <- subset(xenium.obj.15, subset = nCount_Xenium > 0)

path.16 <- "/rsrch3/home/genomic_med/ssatpati/LUNG_Mixed_Histo/20251029__211811__X5K_9MJMHF_Batch2_Lung_10292025/output-XETG00458__0073241__4__20251029__212017/"
#xenium.obj.16 <- LoadXenium(path.16, fov = "fov")
xenium.obj.16 <- ReadXenium(data.dir = path.16,  outs = c("segmentation_method", "matrix", "microns"),  type = "centroids",  mols.qv.threshold = 20)
gene_expr16 <- xenium.obj.16$matrix[["Gene Expression"]]
xenium.obj.16 <- CreateSeuratObject(counts = gene_expr16, assay = "Xenium")
xenium.obj.16 <- subset(xenium.obj.16, subset = nCount_Xenium > 0)

path.17 <- "/rsrch3/home/genomic_med/ssatpati/LUNG_Mixed_Histo/20251114__204637__LUNG_MIXED_X5K_9MJMHF_11142025/output-XETG00458__0067378__1__20251114__204800/"
#xenium.obj.17 <- LoadXenium(path.17, fov = "fov")
xenium.obj.17 <- ReadXenium(data.dir = path.17,  outs = c("segmentation_method", "matrix", "microns"),  type = "centroids",  mols.qv.threshold = 20)
gene_expr17 <- xenium.obj.17$matrix[["Gene Expression"]]
xenium.obj.17 <- CreateSeuratObject(counts = gene_expr17, assay = "Xenium")
xenium.obj.17 <- subset(xenium.obj.17, subset = nCount_Xenium > 0)

path.18 <- "/rsrch3/home/genomic_med/ssatpati/LUNG_Mixed_Histo/20251114__204637__LUNG_MIXED_X5K_9MJMHF_11142025/output-XETG00458__0067378__2__20251114__204800/"
#xenium.obj.18 <- LoadXenium(path.18, fov = "fov")
xenium.obj.18 <- ReadXenium(data.dir = path.18,  outs = c("segmentation_method", "matrix", "microns"),  type = "centroids",  mols.qv.threshold = 20)
gene_expr18 <- xenium.obj.18$matrix[["Gene Expression"]]
xenium.obj.18 <- CreateSeuratObject(counts = gene_expr18, assay = "Xenium")
xenium.obj.18 <- subset(xenium.obj.18, subset = nCount_Xenium > 0)

path.19 <- "/rsrch3/home/genomic_med/ssatpati/LUNG_Mixed_Histo/20251114__204637__LUNG_MIXED_X5K_9MJMHF_11142025/output-XETG00458__0067378__3__20251114__204800/"
#xenium.obj.19 <- LoadXenium(path.19, fov = "fov")
xenium.obj.19 <- ReadXenium(data.dir = path.19,  outs = c("segmentation_method", "matrix", "microns"),  type = "centroids",  mols.qv.threshold = 20)
gene_expr19 <- xenium.obj.19$matrix[["Gene Expression"]]
xenium.obj.19 <- CreateSeuratObject(counts = gene_expr19, assay = "Xenium")
xenium.obj.19 <- subset(xenium.obj.19, subset = nCount_Xenium > 0)

path.20 <- "/rsrch3/home/genomic_med/ssatpati/LUNG_Mixed_Histo/20251114__204637__LUNG_MIXED_X5K_9MJMHF_11142025/output-XETG00458__0067381__1__20251114__204759/"
#xenium.obj.20 <- LoadXenium(path.20, fov = "fov")
xenium.obj.20 <- ReadXenium(data.dir = path.20,  outs = c("segmentation_method", "matrix", "microns"),  type = "centroids",  mols.qv.threshold = 20)
gene_expr20 <- xenium.obj.20$matrix[["Gene Expression"]]
xenium.obj.20 <- CreateSeuratObject(counts = gene_expr20, assay = "Xenium")
xenium.obj.20 <- subset(xenium.obj.20, subset = nCount_Xenium > 0)

path.21 <- "/rsrch3/home/genomic_med/ssatpati/LUNG_Mixed_Histo/20251114__204637__LUNG_MIXED_X5K_9MJMHF_11142025/output-XETG00458__0067381__2__20251114__204800/"
#xenium.obj.21 <- LoadXenium(path.21, fov = "fov")
xenium.obj.21 <- ReadXenium(data.dir = path.21,  outs = c("segmentation_method", "matrix", "microns"),  type = "centroids",  mols.qv.threshold = 20)
gene_expr21 <- xenium.obj.21$matrix[["Gene Expression"]]
xenium.obj.21 <- CreateSeuratObject(counts = gene_expr21, assay = "Xenium")
xenium.obj.21 <- subset(xenium.obj.21, subset = nCount_Xenium > 0)

path.22 <- "/rsrch3/home/genomic_med/ssatpati/LUNG_Mixed_Histo/20251114__204637__LUNG_MIXED_X5K_9MJMHF_11142025/output-XETG00458__0067381__3__20251114__204800/"
#xenium.obj.22 <- LoadXenium(path.22, fov = "fov")
xenium.obj.22 <- ReadXenium(data.dir = path.22,  outs = c("segmentation_method", "matrix", "microns"),  type = "centroids",  mols.qv.threshold = 20)
gene_expr22 <- xenium.obj.22$matrix[["Gene Expression"]]
xenium.obj.22 <- CreateSeuratObject(counts = gene_expr22, assay = "Xenium")
xenium.obj.22 <- subset(xenium.obj.22, subset = nCount_Xenium > 0)

path.23 <- "/rsrch3/home/genomic_med/ssatpati/LUNG_Mixed_Histo/20251204__192825__LUNG_MIXED_X5K_9MJMHF_12042025/output-XETG00458__0073108__1__20251204__192943/"
#xenium.obj.23 <- LoadXenium(path.23, fov = "fov")
xenium.obj.23 <- ReadXenium(data.dir = path.23,  outs = c("segmentation_method", "matrix", "microns"),  type = "centroids",  mols.qv.threshold = 20)
gene_expr23 <- xenium.obj.23$matrix[["Gene Expression"]]
xenium.obj.23 <- CreateSeuratObject(counts = gene_expr23, assay = "Xenium")
xenium.obj.23 <- subset(xenium.obj.23, subset = nCount_Xenium > 0)

path.24 <- "/rsrch3/home/genomic_med/ssatpati/LUNG_Mixed_Histo/20251204__192825__LUNG_MIXED_X5K_9MJMHF_12042025/output-XETG00458__0073108__2__20251204__192943/"
#xenium.obj.24 <- LoadXenium(path.24, fov = "fov")
xenium.obj.24 <- ReadXenium(data.dir = path.24,  outs = c("segmentation_method", "matrix", "microns"),  type = "centroids",  mols.qv.threshold = 20)
gene_expr24 <- xenium.obj.24$matrix[["Gene Expression"]]
xenium.obj.24 <- CreateSeuratObject(counts = gene_expr24, assay = "Xenium")
xenium.obj.24 <- subset(xenium.obj.24, subset = nCount_Xenium > 0)

path.25 <- "/rsrch3/home/genomic_med/ssatpati/LUNG_Mixed_Histo/20251204__192825__LUNG_MIXED_X5K_9MJMHF_12042025/output-XETG00458__0073111__1__20251204__192943/"
#xenium.obj.25 <- LoadXenium(path.25, fov = "fov")
xenium.obj.25 <- ReadXenium(data.dir = path.25,  outs = c("segmentation_method", "matrix", "microns"),  type = "centroids",  mols.qv.threshold = 20)
gene_expr25 <- xenium.obj.25$matrix[["Gene Expression"]]
xenium.obj.25 <- CreateSeuratObject(counts = gene_expr25, assay = "Xenium")
xenium.obj.25 <- subset(xenium.obj.25, subset = nCount_Xenium > 0)

path.26 <- "/rsrch3/home/genomic_med/ssatpati/LUNG_Mixed_Histo/20251204__192825__LUNG_MIXED_X5K_9MJMHF_12042025/output-XETG00458__0073111__2__20251204__192943/"
#xenium.obj.26 <- LoadXenium(path.26, fov = "fov")
xenium.obj.26 <- ReadXenium(data.dir = path.26,  outs = c("segmentation_method", "matrix", "microns"),  type = "centroids",  mols.qv.threshold = 20)
gene_expr26 <- xenium.obj.26$matrix[["Gene Expression"]]
xenium.obj.26 <- CreateSeuratObject(counts = gene_expr26, assay = "Xenium")
xenium.obj.26 <- subset(xenium.obj.26, subset = nCount_Xenium > 0)

path.27 <- "/rsrch3/home/genomic_med/ssatpati/LUNG_Mixed_Histo/20251204__192825__LUNG_MIXED_X5K_9MJMHF_12042025/output-XETG00458__0073111__3__20251204__192943/"
#xenium.obj.27 <- LoadXenium(path.27, fov = "fov")
xenium.obj.27 <- ReadXenium(data.dir = path.27,  outs = c("segmentation_method", "matrix", "microns"),  type = "centroids",  mols.qv.threshold = 10)
gene_expr27 <- xenium.obj.27$matrix[["Gene Expression"]]
xenium.obj.27 <- CreateSeuratObject(counts = gene_expr27, assay = "Xenium")
xenium.obj.27 <- subset(xenium.obj.27, subset = nCount_Xenium > 0)

path.28 <- "/rsrch3/home/genomic_med/ssatpati/LUNG_Mixed_Histo/20251204__192825__LUNG_MIXED_X5K_9MJMHF_12042025/output-XETG00458__0073111__4__20251204__192943/"
#xenium.obj.28 <- LoadXenium(path.28, fov = "fov")
xenium.obj.28 <- ReadXenium(data.dir = path.28,  outs = c("segmentation_method", "matrix", "microns"),  type = "centroids",  mols.qv.threshold = 20)
gene_expr28 <- xenium.obj.28$matrix[["Gene Expression"]]
xenium.obj.28 <- CreateSeuratObject(counts = gene_expr28, assay = "Xenium")
xenium.obj.28 <- subset(xenium.obj.28, subset = nCount_Xenium > 0)

xenium.obj.1$type = "Slide1"
#xenium.obj.2$type = "Slide1"
xenium.obj.3$type = "Slide1"
xenium.obj.4$type = "Slide1"
xenium.obj.5$type = "Slide2"
xenium.obj.6$type = "Slide2"
xenium.obj.7$type = "Slide2"
xenium.obj.8$type = "Slide2"
xenium.obj.9$type = "Slide3"
xenium.obj.10$type = "Slide3"
xenium.obj.11$type = "Slide3"
xenium.obj.12$type = "Slide3"
xenium.obj.13$type = "Slide4"
xenium.obj.14$type = "Slide4"
xenium.obj.15$type = "Slide4"
xenium.obj.16$type = "Slide4"
xenium.obj.17$type = "Slide5"
xenium.obj.18$type = "Slide5"
xenium.obj.19$type = "Slide5"
xenium.obj.20$type = "Slide6"
xenium.obj.21$type = "Slide6"
xenium.obj.22$type = "Slide6"
xenium.obj.23$type = "Slide7"
xenium.obj.24$type = "Slide7"
xenium.obj.25$type = "Slide7"
xenium.obj.26$type = "Slide7"
xenium.obj.27$type = "Slide8"
xenium.obj.28$type = "Slide8"

xenium.obj.1$sampletype = "0067534__1"
#xenium.obj.2$sampletype = "0067534__2"
xenium.obj.3$sampletype = "0067534__3"
xenium.obj.4$sampletype = "0067534__4"
xenium.obj.5$sampletype = "0067555__1"
xenium.obj.6$sampletype = "0067555__2"
xenium.obj.7$sampletype = "0067555__3"
xenium.obj.8$sampletype = "0067555__4"
xenium.obj.9$sampletype = "0073238__1"
xenium.obj.10$sampletype = "0073238__2"
xenium.obj.11$sampletype = "0073238__3"
xenium.obj.12$sampletype = "0073238__4"
xenium.obj.13$sampletype = "0073241__1"
xenium.obj.14$sampletype = "0073241__2"
xenium.obj.15$sampletype = "0073241__3"
xenium.obj.16$sampletype = "0073241__4"
xenium.obj.17$sampletype = "0067378__1"
xenium.obj.18$sampletype = "0067378__2"
xenium.obj.19$sampletype = "0067378__3"
xenium.obj.20$sampletype = "0067381__1"
xenium.obj.21$sampletype = "0067381__2"
xenium.obj.22$sampletype = "0067381__3"
xenium.obj.23$sampletype = "0073108__1"
xenium.obj.24$sampletype = "0073108__2"
xenium.obj.25$sampletype = "0073111__1"
xenium.obj.26$sampletype = "0073111__2"
xenium.obj.27$sampletype = "0073111__3"
xenium.obj.28$sampletype = "0073111__4"

xenium.obj.1$sampleid = "ASC-006-A7"
#xenium.obj.2$sampleid = "ASC-005-A2"
xenium.obj.3$sampleid = "ASC-001-E9"
xenium.obj.4$sampleid = "ASC-010-A4"
xenium.obj.5$sampleid = "ASC-034-B7"
xenium.obj.6$sampleid = "ASC-030-C11"
xenium.obj.7$sampleid = "ASC-031-F6"
xenium.obj.8$sampleid = "ASC-004-F16"
xenium.obj.9$sampleid = "ASC-044-B9"
xenium.obj.10$sampleid = "ASC-037-J4"
xenium.obj.11$sampleid = "ASC-039-A3"
xenium.obj.12$sampleid = "ASC-043-M8"
xenium.obj.13$sampleid = "ASC-036-C5"
xenium.obj.14$sampleid = "ASC-007-I2"
xenium.obj.15$sampleid = "ASC-013-F9"
xenium.obj.16$sampleid = "ASC-008-B11"
xenium.obj.17$sampleid = "ASC-043-M8"
xenium.obj.18$sampleid = "ASC-045-E4"
xenium.obj.19$sampleid = "ASC-053-A4"
xenium.obj.20$sampleid = "ASC-038-D4"
xenium.obj.21$sampleid = "ASC-009-D4"
xenium.obj.22$sampleid = "ASC-040-F17"
xenium.obj.23$sampleid = "ASC-051-E2"
xenium.obj.24$sampleid = "ASC-018-D30"
xenium.obj.25$sampleid = "ASC-048-B8"
xenium.obj.26$sampleid = "ASC-050-B4"
xenium.obj.27$sampleid = "ASC-049-F7"
xenium.obj.28$sampleid = "ASC-049-F7"

#plan("multisession", workers =10)
xenium.obj.1 <- SCTransform(xenium.obj.1, assay = "Xenium", verbose = TRUE, conserve.memory = TRUE)
#xenium.obj.2 <- SCTransform(xenium.obj.2, assay = "Xenium", verbose = FALSE)
xenium.obj.3 <- SCTransform(xenium.obj.3, assay = "Xenium", verbose = TRUE, conserve.memory = TRUE)
xenium.obj.4 <- SCTransform(xenium.obj.4, assay = "Xenium", verbose = TRUE, conserve.memory = TRUE)
xenium.obj.5 <- SCTransform(xenium.obj.5, assay = "Xenium", verbose = TRUE, conserve.memory = TRUE)
xenium.obj.6 <- SCTransform(xenium.obj.6, assay = "Xenium", verbose = TRUE, conserve.memory = TRUE)
xenium.obj.7 <- SCTransform(xenium.obj.7, assay = "Xenium", verbose = TRUE, conserve.memory = TRUE)
xenium.obj.8 <- SCTransform(xenium.obj.8, assay = "Xenium", verbose = TRUE, conserve.memory = TRUE)
xenium.obj.9 <- SCTransform(xenium.obj.9, assay = "Xenium", verbose = TRUE, conserve.memory = TRUE)
xenium.obj.10 <- SCTransform(xenium.obj.10, assay = "Xenium", verbose = TRUE, conserve.memory = TRUE)
xenium.obj.11 <- SCTransform(xenium.obj.11, assay = "Xenium", verbose = TRUE, conserve.memory = TRUE)
xenium.obj.12 <- SCTransform(xenium.obj.12, assay = "Xenium", verbose = TRUE, conserve.memory = TRUE)
xenium.obj.13 <- SCTransform(xenium.obj.13, assay = "Xenium", verbose = TRUE, conserve.memory = TRUE)
xenium.obj.14 <- SCTransform(xenium.obj.14, assay = "Xenium", verbose = TRUE, conserve.memory = TRUE)
xenium.obj.15 <- SCTransform(xenium.obj.15, assay = "Xenium", verbose = TRUE, conserve.memory = TRUE)
xenium.obj.16 <- SCTransform(xenium.obj.16, assay = "Xenium", verbose = TRUE, conserve.memory = TRUE)
xenium.obj.17 <- SCTransform(xenium.obj.17, assay = "Xenium", verbose = TRUE, conserve.memory = TRUE)
xenium.obj.18 <- SCTransform(xenium.obj.18, assay = "Xenium", verbose = TRUE, conserve.memory = TRUE)
xenium.obj.19 <- SCTransform(xenium.obj.19, assay = "Xenium", verbose = TRUE, conserve.memory = TRUE)
xenium.obj.20 <- SCTransform(xenium.obj.20, assay = "Xenium", verbose = TRUE, conserve.memory = TRUE)
xenium.obj.21 <- SCTransform(xenium.obj.21, assay = "Xenium", verbose = TRUE, conserve.memory = TRUE)
xenium.obj.22 <- SCTransform(xenium.obj.22, assay = "Xenium", verbose = TRUE, conserve.memory = TRUE)
xenium.obj.23 <- SCTransform(xenium.obj.23, assay = "Xenium", verbose = TRUE, conserve.memory = TRUE)
xenium.obj.24 <- SCTransform(xenium.obj.24, assay = "Xenium", verbose = TRUE, conserve.memory = TRUE)
xenium.obj.25 <- SCTransform(xenium.obj.25, assay = "Xenium", verbose = TRUE, conserve.memory = TRUE)
xenium.obj.26 <- SCTransform(xenium.obj.26, assay = "Xenium", verbose = TRUE, conserve.memory = TRUE)
xenium.obj.27 <- SCTransform(xenium.obj.27, assay = "Xenium", verbose = TRUE, conserve.memory = TRUE)
xenium.obj.28 <- SCTransform(xenium.obj.28, assay = "Xenium", verbose = TRUE, conserve.memory = TRUE)

#alldata <- merge(xenium.obj.1,c(xenium.obj.2,xenium.obj.3,xenium.obj.4,xenium.obj.5,xenium.obj.6,xenium.obj.7,xenium.obj.8,xenium.obj.9,xenium.obj.10,xenium.obj.11,xenium.obj.12,xenium.obj.13,xenium.obj.14,xenium.obj.15,xenium.obj.16,xenium.obj.17,xenium.obj.18,xenium.obj.19,xenium.obj.20,xenium.obj.21,xenium.obj.22,xenium.obj.23,xenium.obj.24,xenium.obj.25,xenium.obj.26,xenium.obj.27,xenium.obj.28), add.cell.ids=c("0067534__1","0067534__2","0067534__3","0067534__4","0067555__1","0067555__2","0067555__3","0067555__4","0073238__1","0073238__2","0073238__3","0073238__4","0073241__1","0073241__2","0073241__3","0073241__4","0067378__1","0067378__2","0067378__3","0067381__1","0067381__2","0067381__3","0073108__1","0073108__2","0073111__1","0073111__2","0073111__3","0073111__4"))

alldata <- merge(xenium.obj.1,c(xenium.obj.3,xenium.obj.4,xenium.obj.5,xenium.obj.6,xenium.obj.7,xenium.obj.8,xenium.obj.9,xenium.obj.10,xenium.obj.11,xenium.obj.12,xenium.obj.13,xenium.obj.14,xenium.obj.15,xenium.obj.16,xenium.obj.17,xenium.obj.18,xenium.obj.19,xenium.obj.20,xenium.obj.21,xenium.obj.22,xenium.obj.23,xenium.obj.24,xenium.obj.25,xenium.obj.26,xenium.obj.27,xenium.obj.28), add.cell.ids=c("0067534__1","0067534__3","0067534__4","0067555__1","0067555__2","0067555__3","0067555__4","0073238__1","0073238__2","0073238__3","0073238__4","0073241__1","0073241__2","0073241__3","0073241__4","0067378__1","0067378__2","0067378__3","0067381__1","0067381__2","0067381__3","0073108__1","0073108__2","0073111__1","0073111__2","0073111__3","0073111__4"))


# Save the entire R session workspace to an .RData file
save.image("alldata_workspace.RData")

DefaultAssay(alldata) <- "Xenium"

#alldata[["Xenium"]]@counts <- seurat_assay

alldata <- NormalizeData(alldata)
alldata <- ScaleData(alldata)
alldata <- RunPCA(alldata, npcs = 30, features = rownames(alldata))
alldata <- RunUMAP(alldata, dims = 1:30)
alldata <- FindNeighbors(alldata, reduction = "pca", dims = 1:30)
alldata <- FindClusters(alldata, resolution = 0.3)


alldata <- PrepSCTFindMarkers(object = alldata)
alldata <- JoinLayers(alldata)

alldata.markers <- FindAllMarkers(alldata, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
alldata.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
top100_markers <- alldata.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
write.table(top100_markers, file="All_LUNG_100.xls", quote=FALSE, sep="\t", row.names=FALSE)

alldata.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10_markers <- alldata.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.table(top10_markers, file="All_LUNG.markers.xls", quote=FALSE, sep="\t", row.names=FALSE)

names(alldata@meta.data)

table(alldata@meta.data$seurat_clusters)



saveRDS(alldata,'All_Xenium_final.rds')

library(Seurat)
library(ShinyCell)


scConf = createConfig(alldata)
DefaultAssay(alldata) = "SCT"

makeShinyApp(alldata, scConf,gex.assay = "SCT",gene.mapping = TRUE, shiny.dir = "Xenium_LUNG_Shiny", shiny.title = "Dr.Jay_Xenium",shiny.footnotes = "by Rai Lab")

