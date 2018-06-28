rm(list=ls())

#library(cellrangerRkit) # 2.0.0
#library(gtools)
#library(Matrix)
setwd("/mnt/glusterfs/liz86/Wei/SingleCellDE/Code/Function_git")

datawd <- "/mnt/glusterfs/liz86/Wei/scRNASemisupervisedClustering/Data"

# load t cells
load(paste0(datawd, "/cytotoxic_t_matrix.RData"))
dim(data) # 32738 10209
tcells <- data
colnames(tcells) <- paste0(colnames(tcells), "_tcell")
set.seed(20180521)
keep_cell_id <- sample(seq(1, ncol(tcells)), size=100, replace=FALSE)
tcells <- tcells[, keep_cell_id]

# load b cells
load(paste0(datawd, "/b_cells_matrix.RData"))
dim(data) #32738 10085
bcells <- data
set.seed(201805)
keep_cell_id <- sample(seq(1, ncol(bcells)), size=100, replace=FALSE)
bcells <- bcells[, keep_cell_id]
colnames(bcells) <- paste0(colnames(bcells), "_bcell")

# merge 
all(rownames(tcells)==rownames(bcells))
data <- cbind(tcells, bcells)
dim(data) # 32738   200
all_zero_ind <- apply(data, 1, sum)==0
data <- data[!all_zero_ind, ]
dim(data) #9520  200
data[1:5,1:5]
format(object.size(data), units = "Mb", standard = "SI") #15.9 mb

#for (thing in ls()) { message(thing); print(object.size(get(thing)), units='auto') }

cellTypeLables <- c(rep("tcell", ncol(tcells)), rep("bcell", ncol(bcells)))
length(cellTypeLables) # 200

save(data, file="Data/t_b_cells.RData")
save(cellTypeLables, file="Data/cellTypeLables.RData")


