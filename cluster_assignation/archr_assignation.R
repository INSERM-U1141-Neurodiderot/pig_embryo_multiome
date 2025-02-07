library(ArchR)
library(org.Ss.eg.db)
library(SuscrofaTxdb.11.108.july)
library(patchwork)
library(Seurat)
library(dplyr)
library(RColorBrewer)
library(viridis)
library(gprofiler2)

load(file = "/home/adufour/work/rds_storage/omics/archr_all_v7_stemcells.RData")

archrproj <- archrproj_sub

reference <- readRDS("/home/adufour/work/rds_storage/gastrulation/embryo_E11_E15.rds")

# add gene integration matrix
archrproj2 <- addGeneIntegrationMatrix(
    ArchRProj       = archrproj, 
    useMatrix       = "GeneExpressionMatrix",
    matrixName      = "GeneIntegrationMatrix",
    reducedDims     = "Harmony",
    seRNA           = reference,
    addToArrow      = FALSE,
    sampleCellsATAC = 50000,
    sampleCellsRNA  = 50000,
    dimsToUse       = 1:60,
    nGenes          = 4000,
    reduction       = "cca",
    groupRNA        = "Celltypes",
    nameCell        = "predictedCell_Un",
    nameGroup       = "predictedGroup_Un",
    nameScore       = "predictedScore_Un",
    k.anchor        = 80,
    k.filter        = 1600,
    k.score         = 200,
    max.features    = 800,
    n.trees         = 800
)

save.image("/home/adufour/work/rds_storage/omics/assignation_stemcells_arc108.RData")