library(ArchR)
library(org.Ss.eg.db)
library(SuscrofaTxdb.11.108.july)

options(repr.plot.width = 18, repr.plot.height = 17, repr.plot.pointsize = 24)

geneAnnotation <- readRDS("/home/adufour/work/rds_storage/omics/geneannotation.rds")

addArchRThreads(5)
addArchRLocking(locking = TRUE)

excludeChr = c("AEMK02000442.1", "AEMK02000468.1", "AEMK02000631.1", "AEMK02000312.1", "AEMK02000634.1", "AEMK02000418.1", "AEMK02000704.1", "AEMK02000265.1", "AEMK02000673.1", "AEMK02000570.1", "AEMK02000388.1", "AEMK02000628.1", "AEMK02000643.1", "AEMK02000374.1", "AEMK02000676.1", "AEMK02000629.1", "AEMK02000375.1", "AEMK02000551.1", "AEMK02000500.1", "AEMK02000459.1", "AEMK02000582.1", "AEMK02000473.1", "AEMK02000316.1", "AEMK02000555.1", "AEMK02000681.1", "AEMK02000567.1", "AEMK02000237.1", "AEMK02000155.1", "AEMK02000531.1", "AEMK02000491.1", "AEMK02000351.1", "AEMK02000356.1", "AEMK02000319.1", "AEMK02000544.1", "AEMK02000615.1", "AEMK02000346.1", "AEMK02000383.1", "AEMK02000278.1", "AEMK02000621.1", "AEMK02000634.1", "AEMK02000555.1", "AEMK02000418.1", "AEMK02000374.1", "AEMK02000676.1", "AEMK02000237.1", "AEMK02000582.1", "AEMK02000155.1", "AEMK02000615.1", "AEMK02000629.1", "AEMK02000531.1", "AEMK02000704.1", "AEMK02000312.1", "AEMK02000567.1", "AEMK02000351.1", "AEMK02000265.1", "AEMK02000681.1", "AEMK02000316.1", "AEMK02000544.1", "AEMK02000491.1", "AEMK02000459.1", "AEMK02000673.1", "AEMK02000356.1", "AEMK02000570.1", "AEMK02000551.1", "AEMK02000473.1", "AEMK02000319.1", "AEMK02000388.1", "AEMK02000628.1", "AEMK02000643.1", "AEMK02000375.1")
genomeAnnotation <- createGenomeAnnotation(SuscrofaTxdb.11.108.july, standard = FALSE, filter = TRUE, filterChr = excludeChr)

work_dir <- paste0("/home/adufour/work/archr/archr_all_v6")
dir.create(work_dir, recursive = TRUE)
setwd(work_dir)

# Create Arrow file
createArrowFiles(inputFiles  = c("/home/adufour/work/fragencode/workspace/plus4pigs/results/run_cellranger_count/ensembl/arc_j7_108/outs/atac_fragments.tsv.gz", 
                                "/home/adufour/work/fragencode/workspace/plus4pigs/results/run_cellranger_count/ensembl/arc_j9-1_108/outs/atac_fragments.tsv.gz",
                                "/home/adufour/work/fragencode/workspace/plus4pigs/results/run_cellranger_count/ensembl/arc_j9-2_108/outs/atac_fragments.tsv.gz",
                                "/home/adufour/work/fragencode/workspace/plus4pigs/results/run_cellranger_count/ensembl/arc_j11-1_108/outs/atac_fragments.tsv.gz",
                                "/home/adufour/work/fragencode/workspace/plus4pigs/results/run_cellranger_count/ensembl/arc_j11-2_108/outs/atac_fragments.tsv.gz",
                                "/home/adufour/work/fragencode/workspace/plus4pigs/results/run_cellranger_count/ensembl/arc_j11-3_108/outs/atac_fragments.tsv.gz",
                                "/home/adufour/work/fragencode/workspace/plus4pigs/results/run_cellranger_count/ensembl/arc_j11-4_108/outs/atac_fragments.tsv.gz",
                                "/home/adufour/work/fragencode/workspace/plus4pigs/results/run_cellranger_count/ensembl/arc_lw7_108/outs/atac_fragments.tsv.gz",
                                "/home/adufour/work/fragencode/workspace/plus4pigs/results/run_cellranger_count/ensembl/arc_lw9_108/outs/atac_fragments.tsv.gz"), 
                 sampleNames = c("J7", "J9-1", "J9-2", "J11-1", "J11-2", "J11-3", "J11-4", "lw7", "lw9"), 
                 QCDir       = "QualityControl",
                 logFile     = createLogFile(name = "createArrows", 
                                             logDir = "ArchRLogs"),
                 minTSS = 4,
                 minFrags = 1000,
                 maxFrags = 1e+10,
                 GeneScoreMatParams = list(extendUpstream = c(0, 1e+05), extendDownstream = c(0, 1e+05), useTSS = TRUE),
                 excludeChr = "MT",
                 geneAnnotation = geneAnnotation,
                 genomeAnnotation = genomeAnnotation)

archrproj <- ArchRProject(ArrowFiles = c(paste0(work_dir, "/J7.arrow"), paste0(work_dir, "/J9-1.arrow"), paste0(work_dir, "/J9-2.arrow"), paste0(work_dir, "/J11-1.arrow"), paste0(work_dir, "/J11-2.arrow"), paste0(work_dir, "/J11-3.arrow"), paste0(work_dir, "/J11-4.arrow"), paste0(work_dir, "/lw7.arrow"), paste0(work_dir, "/lw9.arrow")),
                geneAnnotation = geneAnnotation,
                genomeAnnotation = genomeAnnotation)

seRNA_1 <- import10xFeatureMatrix(input = "/home/adufour/work/fragencode/workspace/plus4pigs/results/run_cellranger_count/ensembl/arc_j7_108/outs/filtered_feature_bc_matrix.h5", names = "J7", strictMatch = TRUE)
seRNA_2 <- import10xFeatureMatrix(input = "/home/adufour/work/fragencode/workspace/plus4pigs/results/run_cellranger_count/ensembl/arc_j9-1_108/outs/filtered_feature_bc_matrix.h5", names = "J9-1", strictMatch = TRUE)
seRNA_3 <- import10xFeatureMatrix(input = "/home/adufour/work/fragencode/workspace/plus4pigs/results/run_cellranger_count/ensembl/arc_j9-2_108/outs/filtered_feature_bc_matrix.h5", names = "J9-2", strictMatch = TRUE)
seRNA_4 <- import10xFeatureMatrix(input = "/home/adufour/work/fragencode/workspace/plus4pigs/results/run_cellranger_count/ensembl/arc_j11-1_108/outs/filtered_feature_bc_matrix.h5", names = "J11-1", strictMatch = TRUE)
seRNA_5 <- import10xFeatureMatrix(input = "/home/adufour/work/fragencode/workspace/plus4pigs/results/run_cellranger_count/ensembl/arc_j11-2_108/outs/filtered_feature_bc_matrix.h5", names = "J11-2", strictMatch = TRUE)
seRNA_6 <- import10xFeatureMatrix(input = "/home/adufour/work/fragencode/workspace/plus4pigs/results/run_cellranger_count/ensembl/arc_j11-3_108/outs/filtered_feature_bc_matrix.h5", names = "J11-3", strictMatch = TRUE)
seRNA_7 <- import10xFeatureMatrix(input = "/home/adufour/work/fragencode/workspace/plus4pigs/results/run_cellranger_count/ensembl/arc_j11-4_108/outs/filtered_feature_bc_matrix.h5", names = "J11-4", strictMatch = TRUE)
seRNA_8 <- import10xFeatureMatrix(input = "/home/adufour/work/fragencode/workspace/plus4pigs/results/run_cellranger_count/ensembl/arc_lw7_108/outs/filtered_feature_bc_matrix.h5", names = "lw7", strictMatch = TRUE)
seRNA_9 <- import10xFeatureMatrix(input = "/home/adufour/work/fragencode/workspace/plus4pigs/results/run_cellranger_count/ensembl/arc_lw9_108/outs/filtered_feature_bc_matrix.h5", names = "lw9", strictMatch = TRUE)


seRNAcombined<-cbind(assay(seRNA_1), assay(seRNA_2), assay(seRNA_3), assay(seRNA_4), assay(seRNA_5), assay(seRNA_6), assay(seRNA_7), assay(seRNA_8), assay(seRNA_9))
scRNA <- SummarizedExperiment(assays = list(counts = seRNAcombined), rowRanges = rowRanges(seRNA_9))
scRNA@rowRanges@elementMetadata@listData$id <- scRNA@rowRanges@elementMetadata@listData$name

cellsToKeep <- which(getCellNames(archrproj) %in% colnames(scRNA))
archrproj <- subsetArchRProject(ArchRProj = archrproj, cells = getCellNames(archrproj)[cellsToKeep], force = TRUE)

archrproj <- addGeneExpressionMatrix(input = archrproj, verbose = FALSE, seRNA = scRNA)

archrproj <- addCellColData(
  ArchRProj = archrproj,
  data = gsub("-2", "", gsub("-3", "", gsub("-4", "", gsub("-1", "", getCellColData(archrproj)$Sample)))),
  name = "Time",
  cells = getCellNames(archrproj),
  force = FALSE
)

# Filtering doublets
archrproj <- addDoubletScores(archrproj)
archrproj <- filterDoublets(archrproj)

archrproj <- addIterativeLSI(
    ArchRProj = archrproj, 
    clusterParams = list(
      resolution = 2, 
      sampleCells = 30000,
      n.start = 10
    ),
    saveIterations = FALSE,
    useMatrix = "TileMatrix", 
    depthCol = "nFrags",
    name = "LSI_ATAC",
    force = TRUE
)

archrproj <- addIterativeLSI(
    ArchRProj = archrproj, 
    clusterParams = list(
      resolution = 2, 
      sampleCells = 30000,
      n.start = 10
    ),
    saveIterations = FALSE,
    useMatrix = "GeneExpressionMatrix", 
    depthCol = "Gex_nUMI",
    varFeatures = 2500,
    firstSelection = "variable",
    binarize = FALSE,
    name = "LSI_RNA",
    force = TRUE
)

##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                  Combining LSI results for ATAC and RNA                    ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
archrproj <- addCombinedDims(archrproj, reducedDims = c("LSI_ATAC", "LSI_RNA"), name =  "LSI_Combined")

archrproj <- addHarmony(
    ArchRProj = archrproj,
    reducedDims = "LSI_Combined",
    name = "Harmony",
    groupBy = "Sample"
)

##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                            UMAP on the LSI results                         ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
archrproj <- addUMAP(archrproj, reducedDims = "LSI_ATAC", name = "UMAP_ATAC", minDist = 0.8, force = TRUE)
archrproj <- addUMAP(archrproj, reducedDims = "LSI_RNA", name = "UMAP_RNA", minDist = 0.8, force = TRUE)
archrproj <- addUMAP(archrproj, reducedDims = "LSI_Combined", name = "UMAP_Combined", minDist = 0.8, force = TRUE)
archrproj <- addUMAP(archrproj, reducedDims = "Harmony", name = "UMAP_Harmony", nNeighbors = 15, minDist = 0.8, metric = "cosine", force = TRUE)

##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                           Find clusters of cells                           ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##

archrproj <- addClusters(archrproj, reducedDims = "Harmony", name = "Clusters", resolution = 0.1, force = TRUE)

p1 <- plotEmbedding(archrproj, name = "Clusters", embedding = "UMAP_ATAC", size = 2.5, labelAsFactors=F, labelMeans=F)
p2 <- plotEmbedding(archrproj, name = "Clusters", embedding = "UMAP_RNA", size = 2.5, labelAsFactors=F, labelMeans=F)
p3 <- plotEmbedding(archrproj, name = "Clusters", embedding = "UMAP_Combined", size = 2.5, labelAsFactors=F, labelMeans=F)
p4_a <- plotEmbedding(archrproj, name = "Clusters", embedding = "UMAP_Harmony", size = 2.5, labelAsFactors=F, labelMeans=F)

cellsToKeep <- which(getCellNames(archrproj) %in% rownames(getCellColData(archrproj)[which(!getCellColData(archrproj)$Clusters %in% c("C2")),]))
archrproj_sub <- subsetArchRProject(ArchRProj = archrproj, cells = getCellNames(archrproj)[cellsToKeep], force = TRUE)

##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                            UMAP on the LSI results                         ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
archrproj_sub <- addUMAP(archrproj_sub, reducedDims = "LSI_ATAC", name = "UMAP_ATAC", minDist = 0.8, force = TRUE)
archrproj_sub <- addUMAP(archrproj_sub, reducedDims = "LSI_RNA", name = "UMAP_RNA", minDist = 0.8, force = TRUE)
archrproj_sub <- addUMAP(archrproj_sub, reducedDims = "LSI_Combined", name = "UMAP_Combined", minDist = 0.8, force = TRUE)
archrproj_sub <- addUMAP(archrproj_sub, reducedDims = "Harmony", name = "UMAP_Harmony", nNeighbors = 15, minDist = 0.8, metric = "cosine", force = TRUE)

##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                           Find clusters of cells                           ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##

archrproj_sub <- addClusters(archrproj_sub, reducedDims = "Harmony", name = "Clusters", resolution = 0.1, force = TRUE)

p1 <- plotEmbedding(archrproj_sub, name = "Clusters", embedding = "UMAP_ATAC", size = 2.5, labelAsFactors=F, labelMeans=F)
p2 <- plotEmbedding(archrproj_sub, name = "Clusters", embedding = "UMAP_RNA", size = 2.5, labelAsFactors=F, labelMeans=F)
p3 <- plotEmbedding(archrproj_sub, name = "Clusters", embedding = "UMAP_Combined", size = 2.5, labelAsFactors=F, labelMeans=F)
p4_b <- plotEmbedding(archrproj_sub, name = "Clusters", embedding = "UMAP_Harmony", size = 2.5, labelAsFactors=F, labelMeans=F)

##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                       Defining pseudo-bulk replicates                      ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##

archrproj_sub <- addGroupCoverages(ArchRProj = archrproj_sub, groupBy = "Clusters")

pathToMacs2 <- "/usr/local/bioinfo/src/MACS2/venv_MACS-v2.2.7.1/bin/macs2"

archrproj_sub <- addReproduciblePeakSet(
  ArchRProj = archrproj_sub,
  groupBy = "Clusters",
  pathToMacs2 = pathToMacs2,
  genomeSize = 1341049888
)

archrproj_sub <- addPeakMatrix(archrproj_sub)

archrproj_sub <- addPeak2GeneLinks(
    ArchRProj = archrproj_sub,
    reducedDims = "Harmony",
    useMatrix = "GeneExpressionMatrix",
)

pwm_list <- readRDS("/home/adufour/work/rds_storage/omics/pwm_list.rds")

archrproj_sub <- addMotifAnnotations(
  ArchRProj = archrproj_sub,
  name = "Motif",
  motifPWMs = pwm_list,
  cutOff = 5e-05,
  width = 7,
  version = 2,
  force = FALSE,
  logFile = createLogFile("addMotifAnnotations")
)

markerPeaks <- getMarkerFeatures(
  ArchRProj = archrproj_sub, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)", "log10(Gex_nUMI)")
)

saveRDS(markerPeaks, file = paste0("/home/adufour/work/rds_storage/omics/markerPeaks.rds"))

motifsUp <- peakAnnoEnrichment(
    seMarker = markerPeaks,
    ArchRProj = archrproj_sub,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)

saveRDS(motifsUp, file = paste0("/home/adufour/work/rds_storage/omics/motifsUp.rds"))

archrproj_sub <- addBgdPeaks(archrproj_sub)

archrproj_sub <- addDeviationsMatrix(
  ArchRProj = archrproj_sub, 
  peakAnnotation = "Motif"
)

save.image(file = "/home/adufour/work/rds_storage/omics/archr_all_v6.RData")
