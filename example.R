rm(list=ls())

# load B and T cells (downloaded from 10x website)
datawd <- "/net/wong08/home/liz86/SingleCellDE/Code/Function_git/Data/"
load(paste0(datawd, "t_b_cells.RData"))
dim(data) #9466  200
load(paste0(datawd, "cellTypeLables.RData"))
x <- (cellTypeLables=="tcell")*1

exp_folder <- "/net/wong08/home/liz86/SingleCellDE/Code/Function_git/"
code_folder <- "/net/wong08/home/liz86/SingleCellDE/Code/Function_git/code/"

wd <- paste0(exp_folder, "Output/")

if (file.exists(wd)){
    setwd(file.path(wd))
} else {
    dir.create(wd)
    setwd(wd)
}

source(paste0(code_folder, "scRNADEComb.R"))
source(paste0(code_folder, "applyMethods.R"))
TMBfolder <- code_folder

methods <- c("MAST", "MAST_det", "ROTS", "Monocle",
    "DESeq2")
# SCDE failed

p0 <- apply(data, 1, function(x) mean(x!=0))
data <- data[p0 > 0.5, ]
dim(data) #186 200

### apply function #####
ncores <- 20
seed <- 123
sd_above <- NULL
dropout_rate_below <- NULL
log_s_hat <- NULL

ss <- scRNADEComb(count=data, x=x,  
  method=methods, time_name="Comp_time", 
  log_s_hat=NULL, 
  dropout_rate_below=NULL, sd_above=NULL, 
  ncores=ncores, seed=seed, saveDir="res", TMBfolder=TMBfolder)





