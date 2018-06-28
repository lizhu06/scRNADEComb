scRNADEComb <- function(count, x, sample_id=NULL, 
	method=c("MAST", "MAST_det", "ROTS", "Monocle",
		"ZINB", "ZINB_F", "ZINB_R", "SCDE", "DESeq2"), 
	time_name=NULL, log_s_hat = NULL, 
	dropout_rate_below=NULL, sd_above=NULL, 
	ncores=1, seed=123, saveDir=NULL, TMBfolder=NULL, debug=FALSE){

	set.seed(seed)
	options(stringsAsFactors=FALSE)
	# create folder to save results
	mainDir <- getwd()
	if(!is.null(saveDir)){
		if (file.exists(saveDir)){
		    setwd(file.path(mainDir, saveDir))
		} else {
		    dir.create(file.path(mainDir, saveDir))
		    setwd(file.path(mainDir, saveDir))

		}
	}

	# save computing time for each methods
	all_methods <- c("MAST", "MAST_det", "ROTS", "Monocle", 
		"ZINB", "ZINB_F", "ZINB_R", "SCDE", "DESeq2")
	computing_time <- rep(NA, length(all_methods)) 
	names(computing_time) <- all_methods

	# check format/missingness 
	if(!is.matrix(count)){
		count <- data.matrix(count)
	}

	if(is.null(rownames(count))){
		stop("Please use gene names as rownames of count")
	}

	# missing data
	if(any(is.na(count))){stop("Missing data in count.")}
	#if(any(is.na(sample_x))){stop("Missing data in sample x.")}
	if(any(is.na(x))){stop("Missing data in sample x.")}

	# convert x and sample id to factor
	#if(!is.factor(sample_x)){sample_x <- as.factor(sample_x)} # can only take 0/1 (because of MAST)
	if(!is.factor(x)){x <- as.factor(x)} # can only take 0/1 (because of MAST)
	if(!is.null(sample_id)){
		if(!is.factor(sample_id)){sample_id <- as.factor(sample_id)}
	}
	
	#x <- sample_x[sample_id] # FIXME: directly use x as input
	
	if(length(levels(x))<2) {
		stop("Please include at least two groups")
	}

	if(is.null(dropout_rate_below) && is.null(sd_above)){
		message(paste("Total number of genes: ", 
			nrow(count), "\n",
			"Total number of samples: ", ncol(count), "\n",
			"No filtering is applied.\n"))
	}

	# gene filtering
	keep_index1 <- NULL
	keep_index2 <- NULL
	if(!is.null(dropout_rate_below)) {
		dr <- apply(count==0, 1, mean)
		keep_index1 <- dr < dropout_rate_below
	}

	if(!is.null(sd_above)) {
		sd <- apply(count, 1, sd)
		keep_index2 <- sd > sd_above
	}

	if((!is.null(keep_index1)) && is.null(keep_index2)){
		keep_index <- keep_index1
		count <- count[keep_index, ]
		cat(paste("After filtering, ", nrow(count), " genes left. \n", sep=""))
	}else if(is.null(keep_index1) && (!is.null(keep_index2))){
		keep_index <- keep_index2
		count <- count[keep_index, ]
		cat(paste("After filtering, ", nrow(count), " genes left. \n", sep=""))
	}else if((!is.null(keep_index1)) && (!is.null(keep_index2))){
		keep_index <- rep(TRUE, nrow(count))
		count <- count[keep_index, ]
		cat(paste("After filtering, ", nrow(count), " genes left. \n", sep=""))
	}

	# pre-estimate beta_nu_p and beta_s_p for ZINB
	if(length(intersect(c("ZINB", "ZINB_F", "ZINB_R"), method))>0) {
	
		if(!("TMB" %in% rownames(installed.packages())))
		{
			install.packages("TMB", repos='http://cran.us.r-project.org')
			library(TMB)
		}else{
			library(TMB)
		}

		if(is.null(TMBfolder)){stop("Please specify the folder of TMB cpp files")}
		compile(paste0(TMBfolder, "SCMMRan.cpp"))
		compile(paste0(TMBfolder, "SCMM.cpp"))

		# check if experiments is nested, if yes, get sample x. Only matters for ZINB_F
		sample_x <- NULL # input is x, not sample x
		if(!is.null(sample_id)){
			unique_sample_id <- levels(sample_id)
			num_uni_x_per_sample <- sapply(1:length(unique_sample_id), function(s) 
				length(unique(x[sample_id == unique_sample_id[s]])))
			if(all(num_uni_x_per_sample==1)){ 
				nested <- TRUE 
				cat(paste0("This is a nested design. \n"))
			}else{
				nested <- FALSE
				cat(paste0("This is not a nested design. \n"))
			}
			if(nested){ 
				#sample_x <- x[match(unique_sample_id, sample_id)] 
				sample_x <- x[match(levels(sample_id), sample_id)] 
			}else{
				sample_x <- NULL
			}
		}
		
		# calculate log(size factor)
		if(is.null(log_s_hat)){
			total_UMI <- colSums(count)
			# deal with cells with counts as all zero's
			total_UMI[which(total_UMI == 0)] <- min(total_UMI[total_UMI > 0])
			s_hat <- total_UMI/median(total_UMI)
		}
		log_s_hat <- log(s_hat)
		save(log_s_hat, file="log_s_hat.RData")

		if(all(log_s_hat==0)){
			no_log_s_hat <- TRUE
		}else{
			no_log_s_hat <- FALSE
		}

	}

	#### DE ####
	#res <- NULL
	if("MAST" %in% method){
		cat("Start MAST.\n")
		t0 <- proc.time()[3]
		#res <- MAST_f(count, x, ncores)
		res <- tryCatch(MAST_f(count, x, ncores), 
			error=function(e){
				cat(paste0("Error in fitting MAST with following error message: \n", 
					e$message, "\n"))
				return(NA)
			}
		)
		cat("Finish MAST.\n")
		computing_time["MAST"] <- proc.time()[3]-t0
		save(res, file="MAST_res.RData")
	}

	if("MAST_det" %in% method){
		cat("Start MAST_det.\n")
		t0 <- proc.time()[3]
		#res <- MAST_det_f(count, x, sample_id, ncores)
		res <- tryCatch(MAST_det_f(count, x, sample_id, ncores), 
			error=function(e) 
			{
				cat(paste0("Error in fitting MAST_det with following error message: \n", 
					e$message, "\n"))
				return(NA)
			}
		)
		cat("Finish MAST_det.\n")
		computing_time["MAST_det"] <- proc.time()[3]-t0
		save(res, file="MAST_det_res.RData")
	}
	
	if("ROTS" %in% method){
		cat("Start ROTS.\n")
		t0 <- proc.time()[3]
		#res <- ROTS_f(count, x)
		res <- tryCatch(ROTS_f(count, x), 
			error=function(e) {
				cat(paste0("Error in fitting ROTS with following error message: \n", 
					e$message, "\n"))
				return(NA)
			}
		)
		cat("Finish ROTS.\n")
		computing_time["ROTS"] <- proc.time()[3]-t0
		save(res, file="ROTS_res.RData")
	}

	if("Monocle" %in% method){
		cat("Start Monocle.\n")
		t0 <- proc.time()[3]
		#res <- Monocle_f(count, x, ncores)
		res <- tryCatch(Monocle_f(count, x, ncores), 
			error=function(e) {
				cat(paste0("Error in fitting Monocle with following error message: \n", 
					e$message, "\n"))
				return(NA)
			}
		)
		cat("Finish Monocle.\n")
		computing_time["Monocle"] <- proc.time()[3]-t0
		save(res,file="Monocle_res.RData")
	}

	if("ZINB" %in% method){
		cat("Start ZINB.\n")
		t0 <- proc.time()[3]
		#tryCatch(res <- ZINB_f(count, x, ncores), 
		#	error=function(e)cat("Error in fitting ZINB."))
		res <- ZINB_f(count, sample_x, x, sample_id, log_s_hat, no_log_s_hat, 
			subject_effect="none", ncores, TMBfolder, debug)

		cat("Finish ZINB.\n")
		computing_time["ZINB"] <- proc.time()[3]-t0
		save(res, file="ZINB_res.RData")
	}

	if("ZINB_F" %in% method){
		cat("Start ZINB_F.\n")
		t0 <- proc.time()[3]
		#tryCatch(res <- ZINB_f(count, x, ncores), 
		#	error=function(e)cat("Error in fitting ZINB."))
		res <- ZINB_f(count, sample_x, x, sample_id, log_s_hat, no_log_s_hat, 
			subject_effect="fixed", ncores, TMBfolder, debug)
		
		cat("Finish ZINB_F.\n")
		computing_time["ZINB_F"] <- proc.time()[3]-t0
		save(res, file="ZINB_F_res.RData")
	}

	if("ZINB_R" %in% method){
		cat("Start ZINB_R.\n")
		t0 <- proc.time()[3]
		#tryCatch(res <- ZINB_f(count, x, ncores), 
		#	error=function(e)cat("Error in fitting ZINB."))
		res <- ZINB_f(count, sample_x, x, sample_id, log_s_hat, no_log_s_hat, 
			subject_effect="random", ncores, TMBfolder, debug)

 		cat("Finish ZINB_R.\n")
		computing_time["ZINB_R"] <- proc.time()[3]-t0
		save(res, file="ZINB_R_res.RData")
	}

	if("SCDE" %in% method){
		cat("Start SCDE.\n")
		t0 <- proc.time()[3]
		#res <- SCDE_f(count, x, ncores)
		res <- tryCatch(SCDE_f(count, x, ncores), 
			error=function(e) {
				cat(paste0("Error in fitting SCDE with following error message: \n", 
					e$message, "\n"))
				return(NA)
			}
		)
		cat("Finish SCDE.\n")
		computing_time["SCDE"] <- proc.time()[3]-t0
		save(res, file="SCDE_res.RData")
	}

	if( "DESeq2" %in% method){
		cat("Start DESeq2.\n")
		t0 <- proc.time()[3]
		#res <- DESeq2_f(count, x, ncores)
		res <- tryCatch(DESeq2_f(count, x, ncores), 
			error=function(e) {
				cat(paste0("Error in fitting DESeq2 with following error message: \n", 
					e$message, "\n"))
				return(NA)
			}	
		)
		cat("Finish DESeq2.\n")
		computing_time["DESeq2"] <- proc.time()[3]-t0
		save(res,file="DESeq2_res.RData")
	}

	if(is.null(time_name)){
		time_name <- "Computing_time"
	}

	save(computing_time, file=paste(time_name, ".RData", sep=""))
	print(sessionInfo())
	setwd(file.path(mainDir))
	return(computing_time)
}


