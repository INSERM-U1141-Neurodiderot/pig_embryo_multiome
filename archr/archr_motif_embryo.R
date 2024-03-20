library(ArchR)
library(org.Ss.eg.db)
library(SuscrofaTxdb.11.108.july)

options(repr.plot.width = 18, repr.plot.height = 17, repr.plot.pointsize = 24)

load(file = "/home/adufour/work/rds_storage/omics/archr_epiblast_v2.RData")

addArchRThreads(1)
addArchRLocking(locking = TRUE)

pwm_list <- readRDS("/home/adufour/work/rds_storage/omics/pwm_list.rds")

archrproj_sub <- addMotifAnnotations(
  ArchRProj = archrproj_sub,
  name = "Motif",
  motifPWMs = pwm_list,
  cutOff = 5e-05,
  width = 7,
  version = 2,
  force = TRUE,
  logFile = createLogFile("addMotifAnnotations")
)

save.image(file = "/home/adufour/work/rds_storage/omics/archr_epiblast_v2_motif.RData")
