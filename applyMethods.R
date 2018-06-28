SCDE_f <- function(count, x, ncores){
	if(!("scde" %in% rownames(installed.packages())))
	{
	  source("https://bioconductor.org/biocLite.R")
	  biocLite("scde")
	  library(scde)
	}else{
		library(scde)
	}
	count <- apply(count, 2, function(x) {storage.mode(x) <- 
		'integer'; x})
	o.ifm <- scde.error.models(counts=count, 
		groups=x, n.cores=ncores, 
		threshold.segmentation=TRUE, 
		save.crossfit.plots=FALSE, 
		save.model.plots=FALSE, verbose=1, 
		min.size.entries=min(2000, nrow(count) - 1))
	valid.cells <- o.ifm$corr.a > 0
	o.ifm <- o.ifm[valid.cells, ]
	priors <- scde.expression.prior(models=o.ifm, 
		counts=count[,valid.cells], 
		length.out=400, show.plot=FALSE
	)
	x2 <- x[which(valid.cells)]
	names(x2) <- rownames(o.ifm)
	res <- scde.expression.difference(o.ifm,
	  count[,valid.cells], priors, groups = x2,
	  n.randomizations=100, n.cores = ncores,
	  verbose=1)
	res$p.values <- 2*pnorm(abs(res$Z),lower.tail=F)
	return(res)
}

MAST_f <- function(count, x, ncores=1){
	if(!("MAST" %in% rownames(installed.packages())))
	{
		source("https://bioconductor.org/biocLite.R")
		biocLite("MAST")
		library(MAST)
	}else{
		library(MAST)
	}
	data <- log2(count+1)
	cdat <- data.frame(cell_id=colnames(data))
	fdat <- data.frame(feature_id=rownames(data))
	sca <- FromMatrix(data, cdat, fdat)
	
	colData(sca)$condition <- x
	if(ncores > 1){
		options("mc.cores"=ncores)
		zlmCond <- zlm.SingleCellAssay(~ condition, sca, 
			parallel=TRUE)
		options("mc.cores"=1)
	}else{
		zlmCond <- zlm.SingleCellAssay(~ condition, sca, 
			parallel=FALSE)
	}
	summaryCond <- summary(zlmCond, doLRT='condition1') # condition=0 is the ref
	summaryDt <- summaryCond$datatable

	res <- merge(merge(summaryDt[contrast=='condition1' & component=='D',
		.(primerid, `Pr(>Chisq)`, coef)], # discrete p
		summaryDt[contrast=='condition1' & component=='C',
		.(primerid, `Pr(>Chisq)`, coef)],  # continuous
		by='primerid'), summaryDt[contrast=='condition1' & component=='H',
	.(primerid, `Pr(>Chisq)`)],  # Hurdle)
	by='primerid') 
	colnames(res)[2:6] <- c("logit_p", "logit_beta", "normal_p", 
	"normal_beta","hurdle_p")
	res$logit_q <- p.adjust(res$logit_p, method="BH")
	res$normal_q <- p.adjust(res$normal_p, method="BH")
	res$hurdle_q <- p.adjust(res$hurdle_p, method="BH")
	res$geneName <- res$primerid
	res$pval_comp <- res$hurdle_p
	res$qval_comp <- res$hurdle_q
	return(res)
}

MAST_det_f <- function(count, x, sample_id=NULL, ncores=1){
	if(!("MAST" %in% rownames(installed.packages())))
	{
		source("https://bioconductor.org/biocLite.R")
		biocLite("MAST")
		library(MAST)
	}else{
		library(MAST)
	}
	data <- log2(count+1)
	cdat <- data.frame(cell_id=colnames(data))
	fdat <- data.frame(feature_id=rownames(data))
	sca <- FromMatrix(data, cdat, fdat)
	colData(sca)$cngeneson <- scale(colSums(assay(sca)>0))

	colData(sca)$condition <- x
	if(ncores > 1){
		options("mc.cores"=20)
		zlmCond <- zlm.SingleCellAssay(~ condition+cngeneson, sca, 
			parallel=TRUE)
		options("mc.cores"=1)
	}else{
		zlmCond <- zlm.SingleCellAssay(~ condition+cngeneson, sca, 
			parallel=FALSE)
	}
	summaryCond <- summary(zlmCond, doLRT='condition1') 
	summaryDt <- summaryCond$datatable

	res <- merge(merge(summaryDt[contrast=='condition1' & component=='D',
		.(primerid, `Pr(>Chisq)`, coef)], # discrete p
		summaryDt[contrast=='condition1' & component=='C',
		.(primerid, `Pr(>Chisq)`, coef)],  # continuous
		by='primerid'), summaryDt[contrast=='condition1' & component=='H',
	.(primerid, `Pr(>Chisq)`)],  # Hurdle)
	by='primerid') 
	colnames(res)[2:6] <- c("logit_p", "logit_beta", "normal_p", 
	"normal_beta","hurdle_p")
	res$logit_q <- p.adjust(res$logit_p, method="BH")
	res$normal_q <- p.adjust(res$normal_p, method="BH")
	res$hurdle_q <- p.adjust(res$hurdle_p, method="BH")
	res$geneName <- res$primerid
	res$pval_comp <- res$hurdle_p
	res$qval_comp <- res$hurdle_q
	return(res)
}


DESeq2_f <- function(count, x, ncores=1){
	if(!("DESeq2" %in% rownames(installed.packages())))
	{
		source("https://bioconductor.org/biocLite.R")
		biocLite("DESeq2")
		library(DESeq2)
	}else{
		library(DESeq2)
	}
	
	colData <- data.frame(label=x)
	dds <- DESeqDataSetFromMatrix(countData=count,
	       colData=colData, design= ~label)
	dds <- estimateSizeFactors(dds, type="iterate") # befault is not working because each gene has zero
	if(ncores > 1){
		if(!("BiocParallel" %in% rownames(installed.packages())))
		{
			source("https://bioconductor.org/biocLite.R")
			biocLite("BiocParallel")
			library(BiocParallel)
		}else{
			library(BiocParallel)
		}
		register(MulticoreParam(ncores))
		dds <- DESeq(dds, parallel=T)
	}else{
		dds <- DESeq(dds, parallel=FALSE)
	}
	res <- results(dds)
}

ROTS_f <- function(count, x){
	if(!("ROTS" %in% rownames(installed.packages())))
	{
		source("https://bioconductor.org/biocLite.R")
		biocLite("ROTS")
		library(ROTS)
	}else{
		library(ROTS)
	}
	data <- log2(count+1)
	res <- ROTS(data=data, groups=as.numeric(as.character(x)), B=1000, K=500, seed=1234) 
	res2 <- data.frame(logfc=res$logfc, pvalue=res$pvalue, FDR=res$FDR)
	rownames(res2) <- names(res$logfc)
	res2$geneName <- rownames(res2)
	res2$pval_comp <- res2$pval
	res2$qval_comp <- res2$FDR
	return(res2)
}

Monocle_f <- function(count, x, ncores=1){
	if(!("monocle" %in% rownames(installed.packages())))
	{
		source("https://bioconductor.org/biocLite.R")
		biocLite("monocle")
		library(monocle)
	}else{
		library(monocle)
	}
	pd <- new("AnnotatedDataFrame", data=data.frame(x))
	rownames(pd) <- colnames(count)
	fd <- new("AnnotatedDataFrame", data=data.frame(rownames(count)))
	rownames(fd) <- rownames(count)
	cds <- newCellDataSet(count, 
		phenoData=pd, featureData=fd, expressionFamily=negbinomial.size())
	cds <- estimateSizeFactors(cds)
	cds <- estimateDispersions(cds, fitType="local") # may fail in some case
	res <- differentialGeneTest(cds,
	          fullModelFormulaStr="~x",
	          cores=ncores)
	res$geneName <- rownames(res)
	res$pval_comp <- res$pval
	res$qval_comp <- res$qval
	return(res)
}


## ZINB
zinb_single_gene <- function(count, sample_x, x, sample_id, log_s_hat, no_log_s_hat, 
	subject_effect, gene_id, debug){
	if(debug){
		cat(paste0("gene_id=", gene_id, "\n"))
	}

	if(gene_id%%100==0){
		cat(paste0("Start gene ", gene_id, "\n"))
	}
	
	seed <- 2017*gene_id+100
	set.seed(seed)

	y <- count[gene_id,,drop=FALSE]
	#x <- matrix(x, ncol=1)

	if(subject_effect=="none"){
			if(debug){
				fit_zinb <- SCMM.gene(y, x, x, sample_x, offsetMu=log_s_hat, offsetP=NULL,
					sample_id=NULL, ZINB_F=FALSE, ZINB_R=FALSE, debug=debug)
				}else{
					fit_zinb <- tryCatch(SCMM.gene(y, x, x, sample_x, offsetMu=log_s_hat, offsetP=NULL,
						sample_id=NULL, ZINB_F=FALSE, ZINB_R=FALSE, debug=debug),
						error=function(e) cat(paste("error message: ", e$message, "\n")))
				}

	}else if(subject_effect=="fixed"){
			if(debug){
				fit_zinb <- SCMM.gene(y, x, x, sample_x, offsetMu=log_s_hat, offsetP=NULL,
					sample_id=sample_id, ZINB_F=TRUE, ZINB_R=FALSE, debug=debug)
				}else{
					fit_zinb <- tryCatch(SCMM.gene(y, x, x, sample_x, offsetMu=log_s_hat, offsetP=NULL,
						sample_id=sample_id, ZINB_F=TRUE, ZINB_R=FALSE, debug), 
						error=function(e) cat(paste("error message: ", e$message, "\n")))
				}

	}else{
		if(debug){
			fit_zinb <- SCMM.gene(y, x, x, sample_x, offsetMu=log_s_hat, offsetP=NULL,
				sample_id=sample_id, ZINB_F=FALSE, ZINB_R=TRUE, debug=debug)
			}else{
				fit_zinb <- tryCatch(SCMM.gene(y, x, x, sample_x, offsetMu=log_s_hat, offsetP=NULL,
					sample_id=sample_id, ZINB_F=FALSE, ZINB_R=TRUE, debug=debug), 
					error=function(e) cat(paste("error message: ", e$message, "\n")))
			}
		
	}
	colnames1 <- c("nb_beta0", "nb_beta0_p", "nb_beta", "nb_beta_p",
		"nb_overdisp", "sd_mu_random","sd_p_random","logit_beta0", "logit_beta0_p", "logit_beta", "logit_beta_p")
	colnames2 <- paste("nb_beta_sample_id", levels(sample_id)[-1], sep="_")
	colnames3 <- paste("nb_beta_p_sample_id", levels(sample_id)[-1], sep="_")
	colnames4 <- paste("logit_beta_sample_id", levels(sample_id)[-1], sep="_")
	colnames5 <- paste("logit_beta_p_sample_id", levels(sample_id)[-1], sep="_")

	col_names <- c(colnames1, colnames2, colnames3, colnames4, colnames5)
	res <- matrix(NA, 1, length(col_names))
	colnames(res) <- col_names
	if(all(is.na(fit_zinb))){
		res[1, "nb_beta0"] = NA
		res[1, "nb_beta0_p"] = NA
		res[1, "nb_beta"]=NA
		res[1, "nb_beta_p"]=NA
		res[1, "nb_overdisp"]=NA
		res[1, "sd_mu_random"]=NA
		res[1, "sd_p_random"]=NA
		res[1, "logit_beta0"]=NA
		res[1, "logit_beta0_p"]=NA
		res[1,"logit_beta"]=NA
		res[1,"logit_beta_p"]=NA
	}else{
		
		#res_summary <- summary(fit_zinb)
		res[1,"nb_overdisp"] <- ifelse(fit_zinb$theta==0, NA, 1/fit_zinb$theta)  # var=mu+phi*mu^2
		res[1,"sd_mu_random"] <- fit_zinb$sd_mu_random 
		res[1,"sd_p_random"] <- fit_zinb$sd_p_random 
		res[1,"logit_beta0"] <- fit_zinb$logit[1,1]
		res[1,"logit_beta0_p"] <- fit_zinb$logit[1,4]

		res[1,"nb_beta0"] <- fit_zinb$nb[1,1]
		res[1,"nb_beta0_p"] <- fit_zinb$nb[1,4]

		if(subject_effect == "fixed"){
		
			res[1, "nb_beta"] <- fit_zinb$nb["nb_beta", 1]
			res[1, "nb_beta_p"] <- fit_zinb$nb["nb_beta", 4]
			res[1, "logit_beta"] <- fit_zinb$logit["logit_beta", 1]
			res[1, "logit_beta_p"] <- fit_zinb$logit["logit_beta", 4]


			res[1, (length(colnames1)+1) : length(c(colnames1,colnames2))] <- 
				fit_zinb$nb[paste("sample_id", levels(sample_id)[-1], sep=""),1]

			res[1, (length(c(colnames1,colnames2))+1) : length(c(colnames1, colnames2, colnames3))] <- 
				fit_zinb$nb[paste("sample_id", levels(sample_id)[-1], sep=""),4]

			res[1, (length(c(colnames1, colnames2, colnames3))+1) : length(c(colnames1, colnames2, colnames3, colnames4))] <- 
				fit_zinb$logit[paste("sample_id", levels(sample_id)[-1], sep=""),1]

			res[1, (length(c(colnames1, colnames2, colnames3, colnames4))+1) : length(col_names)] <- 
				fit_zinb$logit[paste("sample_id", levels(sample_id)[-1], sep=""),4]	
		}else{
			res[1,"nb_beta"] <- fit_zinb$nb[2,1]
			res[1,"nb_beta_p"] <- fit_zinb$nb[2,4]

			res[1,"logit_beta"] <- fit_zinb$logit[2,1]
			res[1,"logit_beta_p"] <- fit_zinb$logit[2,4]
		}
	}
	return(res)
}

zinb_multi_genes <- function(count, sample_x, x, sample_id, log_s_hat, no_log_s_hat, 
	subject_effect, 
	gene_id_vec, gene_name_vec, TMBfolder, debug){

	# can't link both together
	if(subject_effect=="random"){
		#compile(paste0(TMBfolder, "SCMMRan.cpp"))
		dyn.load(dynlib(paste0(TMBfolder, "SCMMRan")))
	}else{
		#compile(paste0(TMBfolder, "SCMM.cpp"))
		dyn.load(dynlib(paste0(TMBfolder, "SCMM")))
	}

	source(paste0(TMBfolder, "SCMM.TMB.v3.R"))

	res <- t(sapply(1:length(gene_id_vec), function(g) 
		zinb_single_gene(count, sample_x, x, sample_id, log_s_hat, no_log_s_hat, 
		subject_effect, gene_id=gene_id_vec[g], debug)))
	colnames1 <- c("nb_beta0", "nb_beta0_p", "nb_beta", "nb_beta_p",
		"nb_overdisp", "sd_mu_random","sd_p_random","logit_beta0", "logit_beta0_p", "logit_beta", "logit_beta_p")
	colnames2 <- paste("nb_beta_sample_id", levels(sample_id)[-1], sep="_")
	colnames3 <- paste("nb_beta_p_sample_id", levels(sample_id)[-1], sep="_")
	colnames4 <- paste("logit_beta_sample_id", levels(sample_id)[-1], sep="_")
	colnames5 <- paste("logit_beta_p_sample_id", levels(sample_id)[-1], sep="_")
	
	col_names <- c(colnames1, colnames2, colnames3, colnames4, colnames5)
	colnames(res) <- col_names
	rownames(res) <- gene_name_vec
	return(res)
	#save(res, file=paste("zinb_res_part", file_suffix, ".RData" ,sep=""))
}

ZINB_f <- function(count, sample_x, x, sample_id, log_s_hat, no_log_s_hat, 
	subject_effect="none", ncores=1, TMBfolder, debug){
	if(!("snowfall" %in% rownames(installed.packages())))
	{
		install.packages("snowfall", repos='http://cran.us.r-project.org')
		library(snowfall)
	}else{
		library(snowfall)
	}
	if(!("TMB" %in% rownames(installed.packages())))
	{
		install.packages("TMB", repos='http://cran.us.r-project.org')
		library(TMB)
	}else{
		library(TMB)
	}
	if(!("MASS" %in% rownames(installed.packages())))
	{
		install.packages("MASS", repos='http://cran.us.r-project.org')
		library(MASS)
	}else{
		library(MASS)
	}

	# need car to test linear combination
	if(!("car" %in% rownames(installed.packages())))
	{
		install.packages("car", repos='http://cran.us.r-project.org')
		library(car)
	}else{
		library(car)
	}

	fl <- floor(nrow(count)/ncores)  # each core will run 1000 genes
	gene_id_list <- lapply(1:(ncores-1), function(core) 
		seq((core-1)*fl+1, core*fl))
	gene_id_list[[ncores]] <- seq(fl*(ncores-1)+1, nrow(count))
	gene_name_list <- lapply(1:ncores, function(core) 
		rownames(count)[gene_id_list[[core]]])

	if(debug){
		res_return <- lapply(seq(1,ncores), function(core) 
		  zinb_multi_genes(count, sample_x, x, sample_id, log_s_hat, no_log_s_hat, 
		  	subject_effect, 
		  	gene_id_list[[core]], gene_name_list[[core]], TMBfolder, debug))
	}else{
		snowfall::sfInit(parallel=TRUE, type="SOCK", cpus=ncores) 
		snowfall::sfExport("zinb_single_gene","zinb_multi_genes", 
			"count", "sample_x", "x", "sample_id",  "log_s_hat", 
			"no_log_s_hat",
			"subject_effect", "gene_id_list", "gene_name_list", "TMBfolder","debug")
		snowfall::sfLibrary(TMB)
		snowfall::sfLibrary(MASS)
		snowfall::sfLibrary(car)
		res_return <- snowfall::sfClusterApply(seq(1,ncores), function(core) 
		  zinb_multi_genes(count, sample_x, x, sample_id, log_s_hat, no_log_s_hat, 
		  	subject_effect, 
		  	gene_id_list[[core]], gene_name_list[[core]], TMBfolder, debug))
		snowfall::sfStop()
	}

	res <- do.call("rbind", res_return)
	res <- data.frame(res)
	res$nb_beta_q <- p.adjust(res$nb_beta_p, method="BH")
	res$logit_beta_q <- p.adjust(res$logit_beta_p, method="BH")
	min_na <- function(x) {
		if(all(is.na(x))){
			return(NA)
		}else{
			return(min(x, na.rm=TRUE))
		}
	}
	res$minp <- as.numeric(apply(res[,c("nb_beta_p","logit_beta_p")], 1, min_na))
	res$minp_p <- pbeta(res$minp, 1, 2)
	res$minp_q <- p.adjust(res$minp_p, method="BH")
	res$geneName <- rownames(res)
	res$pval_comp <- res$minp_p	
	res$qval_comp <- res$minp_q
	return(res)
}

