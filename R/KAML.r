KAML.version <-
function()
{
cat(paste("#", paste(rep("-", 27), collapse=""), "Welcome to KAML", paste(rep("-", 26), collapse=""), "#", sep=""), "\n")
cat("#    ______ _________ ______  _________                              #\n")
cat("#    ___/ //_/___/   |___/  |/  /___/ /  Kinship Adjusted Mult-Locus #\n")
cat("#    __/ ,<   __/ /| |__/ /|_/ / __/ /                BLUP           #\n")
cat("#    _/ /| |  _/ __| |_/ /  / /  _/ /___         Version: 1.1.0      #\n")
cat("#    /_/ |_|  /_/  |_|/_/  /_/   /_____/", "            _\\\\|//_         #\n")
cat("#  Website: https://github.com/YinLiLin/R-KAML      //^. .^\\\\        #\n")
cat(paste("#", paste(rep("-", 47), collapse=""), "ooO-( (00) )-Ooo", paste(rep("-", 5), collapse=""), "#", sep=""), "\n")
}

KAML.stas.cal <- 
function(
	x1, x2, type=c("cor", "auc", "rmse")
)
{
#--------------------------------------------------------------------------------------------------------#
# x1 is the disease phenotype(0/1)																		 #
# x2 is the predicted GEBV																				 #
# Note that x1 and x2 must in the same order															 #
#--------------------------------------------------------------------------------------------------------#
	switch(
		match.arg(type),
		"cor"={
			res <- cor(x1, x2)
		},
		"auc"={
			if(sum(c(0,1) %in% unique(x1)) != 2) stop("Only two levels(case/1,control/0) are allowed!")
			X <- cbind(x1,x2)
			X <- X[order(X[,2]),]
			N.case <- sum(as.numeric(X[,1]) == 1)
			N.cont <- sum(as.numeric(X[,1]) == 0)
			case.index <- which(as.numeric(X[,1]) == 1)
			case.index.mean <- mean(case.index)
			res <- (1/N.cont)*(case.index.mean-0.5*N.case-0.5)
		},
		"rmse"={
			res <- sqrt(sum((x2 - x1)^2) / length(x2))
		}
	)
	return(res)
}

KAML.Crossprod <- 
function(
	x
){

	r.open <- !inherits(try(Revo.version, silent=TRUE),"try-error")
	if(is.null(x)){
		stop("Please assign the matrix!")
	}
	if(r.open){
		cpd <- crossprod(x)
	}else{
		cpd <- crossprodcpp(x)
	}
	return(cpd)
}

KAML.ginv <- 
function(
	x
){
	r.open <- !inherits(try(Revo.version, silent=TRUE),"try-error")
	if(is.null(x)){
		stop("Please assign the matrix!")
	}
	if(r.open){
		gv <- ginv(x)
	}else{
		gv <- geninv(x)
	}
	return(gv)
}

KAML.Bar <- 
function(
	i, n, type=c("type1", "type2", "type3"), symbol="-", tmp.file=NULL, symbol.head="|", symbol.tail=">" ,fixed.points=TRUE, points=seq(0,100,10), symbol.len=50
)
{
#--------------------------------------------------------------------------------------------------------#
# Object: show the rate of progress for loop 	 														 #
#	 																									 #
# Input:	 																							 #
# i: the current loop number	 																		 #
# n: max loop number	 																				 #
# type: type1 for "for" function, type2 for "mclapply" function											 #
# tmp.file: the opened file of "fifo" function															 #
# symbol: the symbol for the rate of progress	 														 #
# symbol.head: the head for the bar																		 #
# symbol.tail: the tail for the bar																		 #
# fixed.points: whether use the setted points which will be printed	 									 #
# points: the setted points which will be printed	 													 #
# symbol.len: the total length of progress bar		 													 #
#--------------------------------------------------------------------------------------------------------#
	switch(
		match.arg(type), 
		"type1"={
			if(fixed.points){
				point.index <- points
				point.index <- point.index[point.index > floor(100*(i-1)/n)]
				if(floor(100*i/n) %in% point.index){
					if(floor(100*i/n) != max(point.index)){
						print.len <- floor(symbol.len*i/n)
						cat(paste("\r", 
							paste(c(symbol.head, rep("-", print.len), symbol.tail), collapse=""), 
							paste(rep(" ", symbol.len-print.len), collapse=""),
							sprintf("%.2f%%", 100*i/n), sep="")
						)
					}else{
						print.len <- floor(symbol.len*i/n)
						cat(paste("\r", 
							paste(c(symbol.head, rep("-", print.len), symbol.tail), collapse=""), 
							sprintf("%.2f%%", 100*i/n), sep="")
						)
					}
				}
			}else{
				if(i < n){
					print.len <- floor(symbol.len*i/n)
					cat(paste("\r", 
						paste(c(symbol.head, rep("-", print.len), symbol.tail), collapse=""), 
						paste(rep(" ", symbol.len-print.len), collapse=""),
						sprintf("%.2f%%", 100*i/n), sep="")
					)
				}else{
					print.len <- floor(symbol.len*i/n)
					cat(paste("\r", 
						paste(c(symbol.head, rep("-", print.len), symbol.tail), collapse=""), 
						sprintf("%.2f%%", 100*i/n), sep="")
					)
				}
			}
		},
		"type2"={
			if(inherits(parallel:::mcfork(), "masterProcess")) {
				progress <- 0.0
				while(progress < n && !isIncomplete(tmp.file)){
					msg <- readBin(tmp.file, "double")
					progress <- progress + as.numeric(msg)
					print.len <- round(symbol.len * progress / n)
					if(fixed.points){
						if(progress %in% round(points * n / 100)){
							cat(paste("\r", 
							paste(c(symbol.head, rep("-", print.len), symbol.tail), collapse=""), 
							paste(rep(" ", symbol.len-print.len), collapse=""),
							sprintf("%.2f%%", progress * 100 / n), sep=""))
						}
					}else{
						cat(paste("\r", 
						paste(c(symbol.head, rep("-", print.len), symbol.tail), collapse=""), 
						paste(rep(" ", symbol.len-print.len), collapse=""),
						sprintf("%.2f%%", progress * 100 / n), sep=""))
					}
				}
				parallel:::mcexit()
			}
		},
		"type3"={
			progress <- readBin(tmp.file, "double") + 1
			writeBin(progress, tmp.file)
			print.len <- round(symbol.len * progress / n)
			if(fixed.points){
				if(progress %in% round(points * n / 100)){
					cat(paste("\r", 
					paste(c(symbol.head, rep("-", print.len), symbol.tail), collapse=""), 
					paste(rep(" ", symbol.len-print.len), collapse=""),
					sprintf("%.2f%%", progress * 100 / n), sep=""))
				}
			}else{
				cat(paste("\r", 
				paste(c(symbol.head, rep("-", print.len), symbol.tail), collapse=""), 
				paste(rep(" ", symbol.len-print.len), collapse=""),
				sprintf("%.2f%%", progress * 100 / n), sep=""))
			}
		}
	)
}

KAML.Bin <-
function(
	GWAS=NULL, bin.size=1000000, inclosure.size=NULL, bin.num=1, MaxBP=1e10, type=c("exact", "approx")
)
{
#--------------------------------------------------------------------------------------------------------#
# Object: To divide Genome into different bins and pick up the top SNPs									 #
#	 																									 #
# Input:	 																							 #
# GWAS: The results of GWAS, 4 columns are included(SNP names, chr, postion, pvalue)	 				 #
# bin.size: The size for each bin																		 #
# inclosure.size: The total number of SNPs that will be picked up from bins	 							 #
# bin.num: How many top SNPs will be picked up from each bins											 #
# MaxBP: Used to distinguish chromosomes	 															 #
# type: Choose exact or approximate algorithm	 														 #
#--------------------------------------------------------------------------------------------------------#
    if(is.null(GWAS))	return(index=NULL)
	options(warn = -1)
	GWAS <- as.matrix(GWAS)
	max.chr <- max(as.numeric(GWAS[, 2]), na.rm=TRUE)
	if(is.infinite(max.chr))	max.chr <- 0
	map.xy.index <- which(!as.numeric(GWAS[, 2]) %in% c(0 : max.chr))
	if(length(map.xy.index) != 0){
		chr.xy <- unique(GWAS[map.xy.index, 2])
		for(i in 1:length(chr.xy)){
			GWAS[GWAS[, 2] == chr.xy[i], 2] <- max.chr + i
		}
	}
	options(warn = 0)
	GWAS <- matrix(as.numeric(GWAS), nrow(GWAS))
    if(is.null(inclosure.size))	type <- "approx"
	switch(
		match.arg(type), 
		"exact" = {
			bin.size <- bin.size / 2
			GWAS.threshold <- 0.05
			theIndex <- NULL
			Temp <- NULL
			inclosure.size.x <- inclosure.size * bin.num
			filtered.index <- which(GWAS[, 4] <= GWAS.threshold)
			filtered.GWAS <- GWAS[filtered.index, ]
			ordered.index <- order(filtered.GWAS[, 4], decreasing=FALSE)
			ordered.GWAS <- filtered.GWAS[ordered.index, ]
			inclosure.size.loop <- length(seq(1, inclosure.size.x, bin.num))
			for(inc in 1: inclosure.size.loop){
				if(nrow(ordered.GWAS) == 0) break()
				if(bin.num == 1){
					theIndex[inc] <- ordered.index[1]
					index1 <- ordered.GWAS[, 2] == ordered.GWAS[1, 2]
					index2 <- (ordered.GWAS[, 3] >= (ordered.GWAS[1, 3]-bin.size)) & (ordered.GWAS[, 3] <= (ordered.GWAS[1, 3] + bin.size))
					ordered.GWAS <- ordered.GWAS[!(index1 & index2), ]
					ordered.index <- ordered.index[!(index1 & index2)]
				}else{
					#thetop <- as.character(ordered.GWAS[1, 1])
					index1 <- ordered.GWAS[, 2] == ordered.GWAS[1, 2]
					index2 <- (ordered.GWAS[, 3] >= ordered.GWAS[1, 3]-bin.size) & (ordered.GWAS[, 3] <= ordered.GWAS[1, 3] + bin.size)
					ordered.GWAS1 <- ordered.GWAS[(index1 & index2), ]
					ordered.index1 <- ordered.index[(index1 & index2)]
					
					#if the number of SNPs in a bin is less than the bin.num, there are some errors. So:
					iRun <- length(ordered.index1) < bin.num
					while(iRun){
						Temp <- c(Temp,ordered.index1)
						ordered.GWAS <- ordered.GWAS[!(index1 & index2),]
						ordered.index <- ordered.index[!(index1 & index2)]
						index1 <- ordered.GWAS[, 2] == ordered.GWAS[1, 2]
						index2 <- (ordered.GWAS[, 3] >= ordered.GWAS[1, 3]-bin.size) & (ordered.GWAS[, 3] <= ordered.GWAS[1, 3] + bin.size)
						ordered.GWAS1 <- ordered.GWAS[(index1 & index2), ]
						ordered.index1 <- ordered.index[(index1 & index2)]
						iRun <- length(ordered.index1) < bin.num
					}
					theIndex[((inc-1) * bin.num + 1):(inc * bin.num)] <- ordered.index1[1:bin.num]
					index3 <- (ordered.GWAS[, 3] >= (min(ordered.GWAS1[1:bin.num, 3])-bin.size)) & (ordered.GWAS[, 3] <= (max(ordered.GWAS1[1:bin.num, 3]) + bin.size))
					ordered.GWAS <- ordered.GWAS[!(index1 & index3), ]
					ordered.index <- ordered.index[!(index1 & index3)]
				}
			}
			theIndex <- filtered.index[c(theIndex,Temp)]
			theIndex <- theIndex[order(GWAS[theIndex, 4], decreasing=FALSE)]
			theIndex <- theIndex[1:inclosure.size]
		}, 
		"approx" = {
			#Create SNP ID: position + CHR * MaxBP
			ID.GP <- as.numeric(as.vector(GWAS[, 3])) + as.numeric(as.vector(GWAS[, 2])) * MaxBP
			
			#Creat bin ID
			bin.GP <- floor(ID.GP/bin.size)
			binP <- as.matrix(cbind(bin.GP, NA, NA, ID.GP, as.numeric(as.vector(GWAS[, 4]))))
			n <- nrow(binP)
			binP <- binP[order(as.numeric(as.vector(binP[, 5]))), ]  #sort on P alue
			binP <- binP[order(as.numeric(as.vector(binP[, 1]))), ]  #sort on bin
			
			#choose the top number in each bin
			binP[(1 + bin.num):n, 2] <- binP[1:(n-bin.num), 1]
			binP[1:bin.num, 2] <- 0
			binP[, 3] <- binP[, 1]-binP[, 2]
			ID.GP <- binP[binP[, 3]>0, ]

			#Handler of single row and remove SNPs without chromosome and position 
			if(is.null(dim(ID.GP))) ID.GP <- matrix(ID.GP, 1, length(ID.GP))
			ID.GP <- ID.GP[order(as.numeric(as.vector(ID.GP[, 5]))), ]
			if(is.null(dim(ID.GP))) ID.GP <- matrix(ID.GP, 1, length(ID.GP))
			index <- !is.na(ID.GP[, 4])
			ID.GP <- ID.GP[index, 4]

			if(!is.null(inclosure.size)){
				if(!is.na(inclosure.size)){
					avaiable <- min(inclosure.size, length(ID.GP))
					if(avaiable == 0){
						ID.GP <- -1
					}else{
						ID.GP <- ID.GP[1:avaiable]
					}
				}
			}
			
			ID.GI <- as.numeric(as.vector(GWAS[, 3])) + as.numeric(as.vector(GWAS[, 2])) * MaxBP
			theIndex <- which(ID.GI %in% ID.GP)
			theIndex <- theIndex[order(GWAS[theIndex, 4], decreasing=FALSE)]
		}
	)
    return(index=theIndex)
}

KAML.HE <- 
function(
	y, X, K
)
{
#--------------------------------------------------------------------------------------------------------#
# Object: To estimate variance component using GEMMA													 #
# Translated from C++ to R by: Haohao Zhang																 #
# Modified by: Lilin Yin 																				 #
# Note: "rootSolve" package needs to be installed														 #
#																										 #
# Input:	 																							 #
# y: The value vector	 																				 #
# X: The fixed effect(X must contain a column of 1's)													 #
# K: Kinship for all individuals(row and column orders must be the same with y)							 #
#--------------------------------------------------------------------------------------------------------#
	#try(setMKLthreads(1),silent = TRUE)
	n = nrow(K)

	#center the Kinship
	w = rep(1, n)
    Gw = rowSums(K)
    alpha = -1.0 / nrow(K)
    K = K + alpha * (w %*% t(Gw) + Gw %*% t(w))

	CalcVChe <- function(y, X, K){
		n = nrow(K)
		r = n / (n - ncol(X))
		v_traceG = mean(diag(K))
		
		# center and scale K by X
		WtW = KAML.Crossprod(X)
		WtWi = solve(WtW)
		WtWiWt = tcrossprod(WtWi, X)
		GW = K %*% X
		Gtmp = GW %*% WtWiWt
		K = K - Gtmp
		Gtmp = t(Gtmp)
		K = K - Gtmp
		WtGW = crossprod(X, GW)
		GW = crossprod(WtWiWt, WtGW)
		Gtmp = GW %*% WtWiWt
		K = K + Gtmp
		d = mean(diag(K))
		traceG_new = d
		if (d != 0) {
			K = K * 1 / d
		}

		# center y by X, and standardize it to have variance 1 (t(y)%*%y/n=1)
		Wty = crossprod(X, y)
		WtWiWty = solve(WtW, Wty)
		y_scale = y - X %*% WtWiWty
		
		VectorVar <- function(x) {
			m = mean(x)
			m2 = sum(x^2) / length(x)
			return (m2 - m * m)
		}
		
		StandardizeVector <- function(x) {
			m = mean(x)
			v = sum((x - m)^2) / length(x)
			v = v - m * m
			x = (x - m) / sqrt(v)
			return (x)
		}
		
		var_y = VectorVar(y)
		var_y_new = VectorVar(y_scale)

		y_scale = StandardizeVector(y_scale)

		# compute Kry, which is used for confidence interval; also compute q_vec (*n^2)
		Kry = K %*% y_scale - r * y_scale
		q_vec = crossprod(Kry, y_scale)
		
		# compuate yKrKKry, which is used later for confidence interval
		KKry = K %*% Kry
		yKrKKry = crossprod(Kry, KKry)
		d = KAML.Crossprod(Kry)
		yKrKKry = c(yKrKKry, d)

		# compute Sij (*n^2)
		tr = sum(K^2)
		S_mat = tr - r * n

		# compute S^{-1}q
		Si_mat = 1/S_mat

		# compute pve (on the transformed scale)
		pve = Si_mat * q_vec

		# compute q_var (*n^4)
		s = 1
		qvar_mat = yKrKKry[1] * pve
		s = s - pve
		tmp_mat = yKrKKry[2] * s
		qvar_mat = (qvar_mat + tmp_mat) * 2.0

		# compute S^{-1}var_qS^{-1}
		Var_mat = Si_mat * qvar_mat * Si_mat

		# transform pve back to the original scale and save data
		s = 1
		vyNtgN = var_y_new / traceG_new
		vtgvy = v_traceG / var_y
		v_sigma2 = pve * vyNtgN
		v_pve = pve * (vyNtgN) * (vtgvy)
		s = s - pve
		pve_total = pve * (vyNtgN) * (vtgvy)

		d = sqrt(Var_mat)
		v_se_sigma2 = d * vyNtgN
		v_se_pve = d * (vyNtgN) * (vtgvy)
		se_pve_total = Var_mat * (vyNtgN) * (vtgvy) * (vyNtgN) * (vtgvy)
		
		v_sigma2 = c(v_sigma2, s * r * var_y_new)
		v_se_sigma2 = c(v_se_sigma2, sqrt(Var_mat) * r * var_y_new)
		se_pve_total = sqrt(se_pve_total)

		return(v_sigma2)
	}

    # initialize sigma2/log_sigma2
    v_sigma2 = CalcVChe(y=y, X=X, K=K)
	
	log_sigma2 = NULL
    for (i in 1:length(v_sigma2)) {
        if (v_sigma2[i] <= 0){
            log_sigma2[i] = log(0.1)
        } else {
            log_sigma2[i] = log(v_sigma2[i])
        }
    }
	
	LogRL_dev1 <- function(log_sigma2, parms) {
		
		y = parms$y
		X = parms$X
		K = parms$K
		n = nrow(K)
		
		P = K * exp(log_sigma2[1]) + diag(n) * exp(log_sigma2[2])

		# calculate H^{-1}
		P = solve(P)

		# calculate P=H^{-1}-H^{-1}X(X^TH^{-1}X)^{-1}X^TH^{-1}
		HiW = P %*% X
		WtHiW = crossprod(X, HiW)
		WtHiWi = solve(WtHiW)

		WtHiWiWtHi = tcrossprod(WtHiWi, HiW)
		P = P - HiW %*% WtHiWiWtHi

		# calculate Py, KPy, PKPy
		Py = P %*% matrix(y)
		KPy = K %*% Py
		
		# calculate dev1=-0.5*trace(PK_i)+0.5*yPKPy
		c(para1 = (-0.5 * sum(t(P) * K) + 0.5 * crossprod(Py, KPy)) * exp(log_sigma2[1]),
			para2 = (-0.5 * sum(diag(P)) + 0.5 * KAML.Crossprod(Py)) * exp(log_sigma2[2]))
	}

	vg = exp(log_sigma2[1])
	ve = exp(log_sigma2[2])
	delta = ve / vg

	return(list(vg=vg, ve=ve, delta=delta))
}

KAML.Mix <- 
function(
	phe, geno=NULL, CV=NULL, K=NULL, eigen.K=NULL, vc.method="emma", lambda=NULL
)
{
#--------------------------------------------------------------------------------------------------------#
# Object: To solve the mix model and calculate the GEBVs												 #
#																										 #
# Input:	 																							 #
# phe: a value vector(NA is allowed)	 																 #
# geno: genotype in numeric format, pure 0, 1, 2 matrix; m * n, m is marker size, n is reference size	 #
# CV: covariance, design matrix(n * x) for the fixed effects(CV must contain a column of 1's)			 #
# K: Kinship for all individuals(row and column orders must be the same with y)							 #
# eigen.K: a priori calculated eigen for K																 #
# vc.method: the methods for estimating variance component						 						 #
# lambda：ve/vg, ve is the variance of residual, vg is the variance of additive effect					 #
#--------------------------------------------------------------------------------------------------------#
    r.open <- !inherits(try(Revo.version, silent=TRUE),"try-error")
	math.cpu <- Math_cpu_check()
	inf.index <- is.na(phe)
	ref.index <- !is.na(phe)
	refphe <- phe[ref.index]
	y <- as.numeric(as.character(refphe))
	N <- length(phe)
	n <- length(refphe)
	
	if(is.null(CV)){
		X <- matrix(1, N, 1)
	}else{
		X <- CV
	}
	refX <- X
	X <- X[ref.index, , drop=FALSE]
	#there is a error when running in Mcrosoft R open with parallel
	if(vc.method != "ai" | (vc.method == "ai" & math.cpu == 1)){
		if(!is.null(eigen.K)){
			eig <- eigen.K
		}else{
			mkl.cpu <- ifelse((2^(n %/% 1000)) < math.cpu, 2^(n %/% 1000), math.cpu)
			try(setMKLthreads(mkl.cpu), silent=TRUE)
			eig <- eigen((K[ref.index, ref.index]), symmetric=TRUE)
			try(setMKLthreads(math.cpu), silent=TRUE)
		}
	}else{
		eig <-NULL
	}
	if(is.null(lambda)){
		if(vc.method == "brent") {
			reml <- KAML.EIGEN.REML(y, X=X, eigenK=eig)
			lambda <- reml$delta
			LL <- NA
		}
        if(vc.method == "emma") {
			reml <- KAML.EMMA.REML(y, X=X, K=K[ref.index, ref.index])
			lambda <- reml$delta
			LL <- reml$REML
		}
		if(vc.method == "he"){
			reml <- KAML.HE(y, X=X, K=K[ref.index, ref.index])
			lambda <- reml$delta
			LL <- NA
		}
		if(vc.method == "ai"){
			if(math.cpu == 1){
			
				# using HE-AI algorithm
				reml <- KAML.HE(y, X=X, K=K[ref.index, ref.index])
				reml <- KAML.AIEM(y, X=X, eigenK=eig, cpu=math.cpu, start=c(reml$vg, reml$ve), verbose=FALSE)
				lambda <- reml$vc[2] / reml$vc[1]
				LL <- NA
			}else{
			
				# using HE-AI algorithm
				reml <- KAML.HE(y, X=X, K=K[ref.index, ref.index])
				reml <- KAML.AIEM(phe, X=refX, K=K, cpu=math.cpu, start=c(reml$vg, reml$ve), verbose=FALSE)
				beta <- reml$beta
				BLUP.ebv <- reml$u
				LL <- NA
				return(list(beta=beta, ebv=BLUP.ebv, LL=LL, reml=reml))
			}
		}
	}
	U <- eig$vectors * matrix(sqrt(1/(eig$values + lambda)), n, length(eig$values), byrow=TRUE);rm(eig); gc()
	yt <- crossprod(U, y)
    X0t <- crossprod(U, X)
	Xt <- X0t
	X0X0 <- KAML.Crossprod(X0t)
	if(X0X0[1, 1] == "NaN"){
		Xt[which(Xt == "NaN")]=0
		yt[which(yt == "NaN")]=0
		XX=KAML.Crossprod(Xt)
	}
	X0Y <- crossprod(X0t, yt)
	XY <- X0Y
	iX0X0 <- try(solve(X0X0), silent = TRUE)
	if(inherits(iX0X0, "try-error")){
		iX0X0 <- KAML.ginv(X0X0)
	}
	iXX <- iX0X0
	beta <- crossprod(iXX, XY)
	XtimesBetaHat <- X %*% beta
	YminusXtimesBetaHat <- matrix(y)- XtimesBetaHat
	Dt <- crossprod(U, YminusXtimesBetaHat)
    S <- U %*% Dt
    BLUP.ebv <- as.vector(K[, ref.index] %*% S)
	return(list(beta=beta, ebv=BLUP.ebv, S=S, LL=LL, reml=reml))
}

KAML.Kin <- 
function(
	M, weight=NULL, type="center", effect="A", priority=c("speed", "memory"), memo=NULL, SUM=NULL, maxLine=1000
)
{
#--------------------------------------------------------------------------------------------------------#
# Object: To calculate the Vanraden Kinship																 #
#																										 #
# Input:	 																							 #
# M: Genotype in numeric format, pure 0, 1, 2 matrix; m * n, m is marker size, n is population size	 	 #
# (both matrix or big.matrix are allowed for M)															 #
# weight: the weights for all makers, the length of weight is m											 #
# effect: the additive(A), dominant(D), epistasis(AA, AD, DD) effect of SNP								 #
# type: which algorithm will be applied("center", "scale", "vanraden")									 #
# priority: choose the calculation speed or memory														 #
# memo: the names of temporary files																	 #
# SUM: the sum will be used to scale the matrix						 									 #
# maxLine：this parameter can control the memory size when using big.matrix								 #
#--------------------------------------------------------------------------------------------------------#
	if(!effect %in% c("A", "D", "AA", "AD", "DD")) stop("Please select the right effect: 'A', 'D', 'AA', 'AD', 'DD'!")
	if(!type %in% c("scale", "center", "vanraden"))	stop("please select the right kinship algorithm: 'center', 'scale', 'vanraden'!")
	if(!is.null(weight)){
		if(sum(is.na(weight)) != 0)	stop("'NA' is not allowed in weight")
		if(sum(weight < 0) != 0)	stop("Negative value is not allowed in weight")
	}
    if(is.null(dim(M))) M <- t(as.matrix(M))
	K <- 1
	SUMx <- SUM
	if(effect == "AD")	Mx <- M
	go <- TRUE
	while(go){
		if(effect == "AD")	{effect <- "A"; AD <- TRUE}else{AD <- FALSE}
		if(effect %in% c("D", "DD")){
			M[which(M == 2, arr.ind=TRUE)] <- 0
			M <- M + 1
		}
		if(effect == "A" & AD)	effect <- "AD"
		switch(
			match.arg(priority), 
			"speed" = {
				if(!is.matrix(M)) M <- as.matrix(M)
				n <- ncol(M)
				m <- nrow(M)
				Pi <- 0.5 * rowMeans(M)
				if(type=="scale"){
					SD <- as.vector(apply(M, 1, sd))
					SD[SD == 0] <- 1
				}
				M <- M-2 * Pi
				if(type=="scale")	M <- M / SD
				if(is.null(SUM)){
					SUM <- ifelse(type=="vanraden", 2 * sum(Pi * (1-Pi)), m)
				}
				rm("Pi"); gc()
				check.r <- !inherits(try(Revo.version,silent=TRUE),"try-error")
				if(!is.null(weight))	M <- M * sqrt(as.vector(weight))
				
				#crossprodcpp is much faster than crossprod in base R
				K <- K * (KAML.Crossprod(M)/SUM)
			}, 
			"memory" = {
				if(!is.big.matrix(M)) stop("Must be big.matrix.")
				n <- ncol(M)
				m <- nrow(M)
				bac <- paste("Z", memo, ".temp.bin", sep="")
				des <- paste("Z", memo, ".temp.desc", sep="")
				#options(bigmemory.typecast.warning=FALSE)
				Z<-big.matrix(nrow=m, ncol=n, type="double", backingfile=bac, descriptorfile=des, init=0.1)
				Pi <- NULL
				estimate.memory <- function(dat, integer=FALSE, raw=FALSE){
					cells.per.gb <- 2^27  # size of double() resulting in ~1GB of memory use by R 2.15
					dimz <- dat 
					if(length(dimz) == 1) { dimz[2] <- 1 }
					if(length(dimz)>1 & length(dimz)<11 & is.numeric(dimz)) {
						total.size <- as.double(1)
						for(cc in 1:length(dimz)) { total.size <- as.double(total.size * as.double(dimz[cc])) }
						memory.estimate <- as.double(as.double(total.size)/cells.per.gb)
						memory.estimate <- memory.estimate
						if(integer) { memory.estimate <- memory.estimate/2 } else { if(raw) { memory.estimate <- memory.estimate/8 } }
						return(memory.estimate)
					} else {
						# guessing this is a vector
						if(!is.list(dimz) & is.vector(dimz)) { 
						  LL <- length(dimz) 
						  return(estimate.memory(LL, integer=integer, raw=raw))
						} else {
						  warning("tried to estimate memory for object which is neither a vector, pair of dimension sizes or a dataframe/matrix") 
						}
					}
				}
				if((Sys.info()[['sysname']]) == 'Windows'){
					max.gb <- memory.limit()/1000
				}else{
					max.gb <- Inf
				}	
				maxLine.gb <- estimate.memory(c(maxLine, n))
				if(maxLine.gb > max.gb) stop("Memory limited! Please reset the 'maxLine'")
				loop.index <- seq(0, m, maxLine)[-1]
				if(max(loop.index) < m) loop.index <- c(loop.index, m)
				loop.len <- length(loop.index)
				cat(" Z assignment...\n")
				for(cc in 1:loop.len){
					if(loop.len == 1){
						c1 <- 1
					}else{
						c1 <- ifelse(cc == loop.len, (loop.index[cc-1]) + 1, loop.index[cc]-maxLine + 1)
					}
					c2 <- loop.index[cc]
					means <-rowMeans(M[c1:c2, 1:n])
					if(!is.null(weight)){
						if(type=="scale"){
							sds <- apply(M[c1:c2, 1:n], 1, sd)
							sds[sds == 0] <- 1
							Z[c1:c2, 1:n] <- (M[c1:c2, 1:n]-means) * sqrt(weight[c1:c2]) / sds
						}else{
							Z[c1:c2, 1:n] <- (M[c1:c2, 1:n]-means) * sqrt(weight[c1:c2])
						}
					}else{
						if(type=="scale"){
							sds <- apply(M[c1:c2, 1:n], 1, sd)
							sds[sds == 0] <- 1
							Z[c1:c2, 1:n] <- (M[c1:c2, 1:n]-means) / sds
						}else{
							Z[c1:c2, 1:n] <- M[c1:c2, 1:n]-means
						}
					}
					Pi <- c(Pi, 0.5 * means);gc()
				}
				cat(" Assignment done!\n")
				if(is.null(SUM)){
					SUM <- ifelse(type=="vanraden", 2 * sum(Pi * (1-Pi)), m)
				}
				fl.suc <- flush(Z)
				if(!fl.suc){ stop("flush failed\n") } 
				RR <- describe(Z); rm(list=c("Z", "Pi")); gc()
				Z <- attach.big.matrix(RR)
				K <- K * (big.crossprod(Z)/SUM)
				rm(Z); gc()
				unlink(c(paste("Z", memo, ".temp.bin", sep=""), paste("Z", memo, ".temp.desc", sep="")), recursive = TRUE)
			}
		)
		if(effect == "AD"){
			effect <- "D"; go <- TRUE;	SUM <- SUMx; M <- Mx
		}else{
			go <- FALSE
		}
	}
	if(effect %in% c("AA", "DD"))	K <- K * K
	return(K)
}

KAML.KinNew <- 
function(
	M, weight=NULL, SUM=NULL, scale=FALSE, priority=c("speed", "memory"), verbose=FALSE, threads=1
){
	if(!is.big.matrix(M))
		stop("genotype must be in 'big.matrix' format.")
	if(!is.null(weight)){
		if(sum(is.na(weight)) != 0)	stop("'NA' is not allowed in weight")
		if(sum(weight < 0) != 0)	stop("Negative value is not allowed in weight")
	}
	priority <- match.arg(priority)
	r.open <- !inherits(try(Revo.version, silent=TRUE),"try-error")
	if(priority == "speed"){
		k <- kin_cal_s(M@address, SUM=SUM, scale=scale, wt=weight, mkl=r.open, threads=threads, verbose=verbose)
	}else{
		k <- kin_cal_m(M@address, SUM=SUM, scale=scale, wt=weight, threads=threads, verbose=verbose)
	}
}

KAML.EMMA.REML <- 
function(
	y, X, K, ngrids=100, llim=-10, ulim=10, esp=1e-10, cpu=1
)
{
#--------------------------------------------------------------------------------------------------------#
# Object: To estimate variance components using optimized EMMA											 #
#																										 #
# Input:	 																							 #
# y: The value vector	 																				 #
# X: The fixed effect(X must contain a column of 1's)													 #
# K: Kinship for all individuals(row and column orders must be the same with y)							 #
# esp: the desired accuracy																				 #
# cpu: the CPU number for calculation																	 #
#--------------------------------------------------------------------------------------------------------#
	R.ver <- Sys.info()[['sysname']]
	wind <- R.ver == 'Windows'
	r.open <- !inherits(try(Revo.version,silent=TRUE),"try-error")
	
	if(!is.numeric(y))	y <- as.numeric(as.character(y))
	emma.delta.REML.LL.wo.Z <- function(logdelta, lambda, etas) 
	{
		nq <- length(etas)
		delta <-  exp(logdelta)
		return( 0.5 * (nq * (log(nq/(2 * pi))-1-log(sum(etas * etas/(lambda + delta))))-sum(log(lambda + delta))) )
	}
	emma.eigen.R.wo.Z=function(K, X) {
		n <- nrow(X)
		q <- ncol(X)
		cXX <- KAML.Crossprod(X)
		iXX <- try(solve(cXX), silent = TRUE)
		if(inherits(iXX, "try-error")){
			iXX <- KAML.ginv(cXX)
		}
		multiply <- function(X, Y){
			myf1 <- function(x){
				apply(x * Y, 2, sum)
			}
			apply(X, 1, myf1)
		}
		SS1 <- X %*% iXX
        SS2 <- tcrossprod(SS1, X)
		S <- diag(n)-SS2
		math.cpu <- try(getMKLthreads(), silent=TRUE)
		mkl.cpu <- ifelse((2^(n %/% 1000)) < math.cpu, 2^(n %/% 1000), math.cpu)
		try(setMKLthreads(mkl.cpu), silent=TRUE)
		eig <- eigen(S %*% (K + diag(1, n)) %*% S, symmetric=TRUE)#S4 error here
		try(setMKLthreads(math.cpu), silent=TRUE)
		stopifnot(!is.complex(eig$values))
		return(list(values=eig$values[1:(n-q)]-1, vectors=eig$vectors[, 1:(n-q)]))
	}
	eig.R=emma.eigen.R.wo.Z(K=K, X=X)
	n <- length(y)
	t <- nrow(K)
	q <- ncol(X)
	#  stopifnot(nrow(K) == t)
	stopifnot(ncol(K) == t)
	stopifnot(nrow(X) == n)
	cXX <- KAML.Crossprod(X)
	if( det(cXX == 0 )){
		warning("X is singular")
		return (list(REML=0, delta=0, ve=0, vg=0))
	}
	etas <- crossprod(eig.R$vectors, y)
	logdelta <- (0:ngrids)/ngrids * (ulim-llim) + llim
	m <- length(logdelta)
	delta <- exp(logdelta)
	Lambdas <- matrix(eig.R$values, n-q, m) + matrix(delta, n-q, m, byrow=TRUE)
	Etasq <- matrix(etas * etas, n-q, m)
	LL <- 0.5 * ((n-q) * (log((n-q)/(2 * pi))-1-log(colSums(Etasq/Lambdas)))-colSums(log(Lambdas)))
	dLL <- 0.5 * delta * ((n-q) * colSums(Etasq/(Lambdas * Lambdas))/colSums(Etasq/Lambdas)-colSums(1/Lambdas))
	optlogdelta <- vector(length=0)
	optLL <- vector(length=0)
	if( dLL[1] < esp ) {
		optlogdelta <- append(optlogdelta, llim)
		optLL <- append(optLL, emma.delta.REML.LL.wo.Z(llim, eig.R$values, etas))
	}
	if( dLL[m-1] > 0-esp ) {
		optlogdelta <- append(optlogdelta, ulim)
		optLL <- append(optLL, emma.delta.REML.LL.wo.Z(ulim, eig.R$values, etas))
	}	
	emma.delta.REML.dLL.wo.Z <- function(logdelta, lambda, etas) {
				nq <- length(etas)
				delta <- exp(logdelta)
				etasq <- etas * etas
				ldelta <- lambda + delta
				return( 0.5 * (nq * sum(etasq/(ldelta * ldelta))/sum(etasq/ldelta)-sum(1/ldelta)) )
	}
	emma.delta.REML.LL.wo.Z <- function(logdelta, lambda, etas) {
				nq <- length(etas)
				delta <-  exp(logdelta)
				return( 0.5 * (nq * (log(nq/(2 * pi))-1-log(sum(etas * etas/(lambda + delta))))-sum(log(lambda + delta))) )
			}
	emma.multi <- function(i){
		if( ( dLL[i] * dLL[i + 1] < 0 ) && ( dLL[i] > 0 ) && ( dLL[i + 1] < 0 ) ) {
			r <- uniroot(emma.delta.REML.dLL.wo.Z, lower=logdelta[i], upper=logdelta[i + 1], lambda=eig.R$values, etas=etas)
			optlogdelta.mult <- r$root
			optLL.mult <- emma.delta.REML.LL.wo.Z(r$root, eig.R$values, etas)
			return(list(optlogdelta.mult=optlogdelta.mult,optLL.mult=optLL.mult))
		}else{
			return(NULL)
		}
	}
	if(cpu > 1){
		if(wind){
			cl <- makeCluster(cpu)
			registerDoParallel(cl)
			clusterExport(cl, varlist=c("emma.delta.REML.dLL.wo.Z", "emma.delta.REML.LL.wo.Z"),envir=environment())
			multi.res <- do.call(rbind, foreach(x=1:(m-1)) %dopar% emma.multi(x))
			stopCluster(cl)
		}else{
			multi.res <- do.call(rbind, parallel::mclapply(1:(m-1), emma.multi, mc.cores=cpu))
		}
	}else{
		multi.res <- do.call(rbind, lapply(1:(m-1), emma.multi))	
	}
	if(!is.null(multi.res)){
		optlogdelta <- c(optlogdelta, unlist(multi.res[,1]))
		optLL <- c(optLL, unlist(multi.res[,2]))
	}
	maxdelta <- exp(optlogdelta[which.max(optLL)])
	#handler of grids with NaN log
	replaceNaN<-function(LL) {
		index=(LL == "NaN")
		if(length(index)>0) theMin=min(LL[!index])
		if(length(index)<1) theMin="NaN"
			LL[index]=theMin
			return(LL)    
	}
	optLL=replaceNaN(optLL)   
	maxLL <- max(optLL)
    maxva <- sum(etas * etas/(eig.R$values + maxdelta))/(n-q)    
	maxve <- maxva * maxdelta
	return (list(REML=maxLL, delta=maxdelta, ve=maxve, vg=maxva))
}

KAML.Impute <- 
function(
	bigm, threads=0
){
	if(is.character(bigm)){
		bigm=attach.big.matrix(paste0(bigm,".geno.desc"))
	}
	if(hasNA(bigm@address)){
		cat(" Imputing missings...\n")
		impute_marker(bigm@address, threads=threads, verbose=FALSE)
	}
}

KAML.Data <-
function(
	hfile=NULL, vfile=NULL, numfile=NULL, mapfile=NULL, bfile=NULL, out=NULL, sep="\t", SNP.impute=c("Left", "Middle", "Right"), maxLine=10000, priority="memory"
)
{
#--------------------------------------------------------------------------------------------------------#
# Object: To prepare data for KAML package																 #
# For binary file:																						 #
# source("https://bioconductor.org/biocLite.R")															 #
# biocLite("snpStats")																					 #
# 	 																									 #
# Input:	 																							 #
# hfile: Genotype in hapmap format	 																	 #
# numfile: Genotype in numeric format; pure 0, 1, 2 matrix; m * n, m is marker size, n is sample size	 #
# mapfile: the map for numfile; three columns: snp, chrom, position										 #
# (Note: hfile and file Num can not input at the same time; numfile and mapfile are both needed)		 #
# bfile: plink binary file, including three files: bfile.fam, bfile.bed, bfile.bim						 #
# out: Name of output files																				 #
# sep: seperator for hapmap, numeric, map, and phenotype												 #
# SNP.impute: "Left", "Middle", "Right"								 									 #
# maxLine：this parameter can control the memory size when using big.matrix								 #
# priority: choose take 'speed' or 'memory' into consideration											 #
#--------------------------------------------------------------------------------------------------------#	
	
	KAML.version()
	cat(" Preparing data for KAML...\n")
	#if(is.null(hfile)&is.null(numfile))
	#stop("Hapmap or Numeric genotype is needed.")
	if(!is.null(hfile)&!is.null(numfile))
	stop("Only one of Hapmap or Numeric genotype is needed!")
	if((!is.null(numfile) & is.null(mapfile)) | (is.null(numfile) & !is.null(mapfile)))
	stop("Both Map and Numeric genotype are needed!")
	if(is.null(out)) out="KAML"
	SNP.impute <- match.arg(SNP.impute)

    if(!is.null(vfile)){
		VCF.Numeralization <- function(
			x, impute="Middle"
		){
			x[x=="0/0"] = 0
			x[x=="0/1"] = 1
			x[x=="1/0"] = 1
			x[x=="1/1"] = 2
			x[x=="./1"] = "N"
			x[x=="1/."] = "N"
			x[x=="./0"] = "N"
			x[x=="0/."] = "N"
			x[x=="./."] = "N"

			#Imputation for missing value
			if(impute=="Middle"){
				x[x=="N"] = 1
			}
			
			if(impute=="Left"){
				n=length(x)
				lev=levels(as.factor(x))
				lev=setdiff(lev,"N")
				len=length(lev)
				if(len==0){
					minA = 1
				}else if(len==1){
					minA = lev[1]
				}else if(len==2){
					len.1 <- length(x[x==lev[1]])
					len.2 <- length(x[x==lev[2]])
					if(len.1<len.2){
						minA = lev[1]
					}else{
						minA = lev[2]
					}
				}else if(len==3){
					len.1 <- length(x[x==lev[1]])
					len.2 <- length(x[x==lev[2]])
					len.3 <- length(x[x==lev[3]])

					min.all = min(len.1, len.2, len.3)
					if(len.1 == min.all){
						minA = lev[1]
					}
					if(len.2 == min.all){
						minA = lev[2]
					}
					if(len.3 == min.all){
						minA = lev[3]
					}
				}
				x[x=="N"] = minA
			}
		
			if(impute=="Right"){
				n=length(x)
				lev=levels(as.factor(x))
				lev=setdiff(lev,"N")
				len=length(lev)
				if(len==0){
					maxA = 1
				}else if(len==1){
					maxA = lev[1]
				}else if(len==2){
					len.1 <- length(x[x==lev[1]])
					len.2 <- length(x[x==lev[2]])
					if(len.1<len.2){
						maxA = lev[2]
					}else{
						maxA = lev[1]
					}
				}else if(len==3){
					len.1 <- length(x[x==lev[1]])
					len.2 <- length(x[x==lev[2]])
					len.3 <- length(x[x==lev[3]])

					max.all = max(len.1, len.2, len.3)

					if(len.1 == max.all){
						maxA = lev[1]
					}
					if(len.2 == max.all){
						maxA = lev[2]
					}
					if(len.3 == max.all){
						maxA = lev[3]
					}
				}
				x[x=="N"] = maxA
			}
			return(x)
		}
    
		dofile <- TRUE
		i <- 0
		fileVCFCon <- file(description=vfile[1], open="r")
		while(dofile){
			i <- i + 1
			tt <- readLines(fileVCFCon, n=1)
			char12 <- substr(tt, 1, 2)
			if(char12 != "##"){
				dofile <- FALSE
				vcf.jump <- i - 1
			}
		}
		cat(" Skip row number:", vcf.jump, "\n")
		close.connection(fileVCFCon)

        #get the first vcf
        fileVCFCon <- file(description=vfile[1], open="r")     
		jj <- readLines(fileVCFCon, n=vcf.jump)
        tt <- readLines(fileVCFCon, n=1, skipNul=1)
        close.connection(fileVCFCon)
		
        #tt2 <- unlist(strsplit(tt, sep))
        tt3 <- unlist(strsplit(tt, sep))
        tt2 <- unlist(strsplit(tt3, "-9_"))
        taxa.g <- as.vector(tt2[-c(1:9)])
        taxa.g <- taxa.g[taxa.g!=""]
        taxa=taxa.g
		
        ns = length(taxa)  #Number of individuals
		cat(" Number of individuals:", ns, "\n")
        nFile=length(vfile)
        
        #Iteration among file
        cat(" Numericilization...\n")
        #if(priority == "memory"){
            for (theFile in 1:nFile){
			
                #Open VCF files
                fileVCFCon <- file(description=vfile[theFile], open="r")

                #handler of first file
                if(theFile == 1){
                    #Open GD and GM file
                    fileNumCon <- file(description=paste(out, ".Numeric.txt", sep=""), open="w")
                    fileMapCon <- file(description=paste(out, ".map", sep=""), open="w")
                }
                
                #Initialization for iteration within file
                inFile=TRUE
                i=0
                jj <- readLines(fileVCFCon, n=vcf.jump+1)
                #Iteration within file
                while(inFile){
                    i =i + 1
                    if(i %% 1000 == 0) cat(paste(" Number of Markers Written into File: ", theFile, ": ", i, sep=""), "\n")

                    tt <- readLines(fileVCFCon, n=1)
                    tt2 <- unlist(strsplit(x=tt, split=sep, fixed=TRUE))
                    
                    #Identify end of file
                    if(is.null(tt2[1])){
						cat(paste(" Number of Markers Written into File: ", theFile, ": ", i-1, sep=""), "\n")
						inFile=FALSE
					}
                    
                    if(inFile){
                        #GM
                        rs=tt2[3]
                        chrom=tt2[1]
                        pos=tt2[2]
                        writeLines(as.character(rs), fileMapCon, sep=sep)
                        writeLines(as.character(chrom), fileMapCon, sep=sep)
                        writeLines(as.character(pos), fileMapCon, sep="\n")
 
                        #GD
                        GD <- VCF.Numeralization(x=tt2[-c(1:9)], impute=SNP.impute)
                        writeLines(as.character(GD[1:(ns-1)]), fileNumCon, sep=sep)
                        writeLines(as.character(GD[ns]), fileNumCon, sep="\n")
                    }#enf of inFIle
                } #end iteration within file
                
                #Close VCF file
				close.connection(fileVCFCon)
            } #end iteration among files
            
            #Close GD and GM file
            close.connection(fileNumCon)
            close.connection(fileMapCon)
        #}
        cat(" Preparation for NUMERIC data is done!\n")
    }
	
	if(!is.null(bfile)){
		map <- read.table(paste(bfile, ".bim", sep=""), head=FALSE, stringsAsFactors=FALSE)
		map <- map[, c(2, 1, 4)]
		cat(" Reading binary files...\n")
		write.table(map, paste(out, ".map", sep=""), row.names=FALSE, col.names=FALSE, sep=sep, quote=FALSE)
		fam <- read.table(paste(bfile, ".fam", sep=""), head=FALSE, stringsAsFactors=FALSE)
		bck <- paste(out, ".geno.bin", sep="")
		dsc <- paste(out, ".geno.desc", sep="")
		nmarkers <- nrow(map)
		ns <- nrow(fam)
		cat(paste(" Number of Ind: ", ns, "; Number of markers: ", nmarkers, ";\n", sep=""))
		myGeno.backed<-big.matrix(nmarkers, ns, type="char",
			backingfile=bck, descriptorfile=dsc)
		cat(" Output BIG genotype...\n")
		TransData_c(bfile = bfile, pBigMat = myGeno.backed@address, maxLine = maxLine, threads = 0)
		KAML.Impute(myGeno.backed)
		# inGENOFile=TRUE
		# i <- 0
		# #printN <- unlist(strsplit(x=as.character(nmarkers), split="", fixed=TRUE))
		# #printIndex <- seq(0, (as.numeric(printN[1]) + 1) * 10^(length(printN)), 1000)[-1]
		# Num.fun <- function(x){
		# 	x <- data.matrix(as.data.frame(x))
		# 	x[x==0]=2
		# 	return(x)
		# }
		# while(inGENOFile){
		# 	i <- i + maxLine
		# 	if(i >= nmarkers){
		# 		xx <- nmarkers
		# 		inGENOFile <- FALSE
		# 	}else{
		# 		xx <- i
		# 	}
		# 	#if(sum(i >= printIndex )!=0){
		# 	#	printIndex <- printIndex[printIndex > i]
		# 	#	print(paste("Number of Markers Written into BIG File: ", xx, sep=""))
		# 	#}
		# 	if(i >= nmarkers){
		# 		myGeno.backed [(i-maxLine + 1):nmarkers, ] <- -1 * apply(geno[, (i-maxLine + 1):nmarkers], 1, Num.fun) + 3
		# 	}else{
		# 		myGeno.backed [(i-maxLine + 1):i, ] <- -1 * apply(geno[, (i-maxLine + 1):i], 1, Num.fun) + 3
		# 	}
		# 	KAML.Bar(i=xx, n=nmarkers, fixed.points=FALSE)
		# }
		geno.flush <- flush(myGeno.backed)
		if(!geno.flush){
			stop("flush failed")
		}else{
			cat(" Preparation for BIG data is done!\n")
		}
		rm(myGeno.backed)
	}

	if(!is.null(hfile)){

		nFile=length(hfile)

		#Iteration among file
		cat(" Output NUMERIC genotype...\n")
		if(priority == "memory"){
			for (theFile in 1:nFile){

				#Open HMP files  
				fileHMPCon<-file(description=hfile[theFile], open="r")
				  
				#Get HMP hearder
				tt<-readLines(fileHMPCon, n=1)
				tt2<-unlist(strsplit(x=tt, split=sep, fixed=TRUE))
				taxa.g<-as.vector(tt2[-c(1:11)] )
				ns<-length(taxa.g)  #Number of individuals

				#handler of first file
				if(theFile == 1){
					#Open GD and GM file
					fileNumCon<-file(description=paste(out, ".Numeric.txt", sep=""), open="w")
					fileMapCon<-file(description=paste(out, ".map", sep=""), open="w")
					#GM header
					# writeLines("SNP", fileMapCon, sep=sep)
					# writeLines("Chrom", fileMapCon, sep=sep)
					# writeLines("BP", fileMapCon, sep="\n")  
				}
		  
				#Initialization for iteration within file
				inFile=TRUE
				i=0
				#Iteration within file
				while(inFile){
					i=i + 1
					if(i %% 1000 == 0)cat(paste(" Number of Markers finished for theFile ", theFile, ": ", i, sep=""), "\n")
					tt<-readLines(fileHMPCon, n=1) 
					tt2<-unlist(strsplit(x=tt, split=sep, fixed=TRUE))
				
					#Identify end of file
					if(is.null(tt2[1]))inFile=FALSE
					
					if(inFile){
						#GM
						rs=tt2[1]
						chrom=tt2[3]
						pos=tt2[4]
						writeLines(as.character(rs), fileMapCon, sep=sep)
						writeLines(as.character(chrom), fileMapCon, sep=sep)
						writeLines(as.character(pos), fileMapCon, sep="\n")
					  
						#GD
						GD = KAML.Num(x=tt2[-c(1:11)], impute=SNP.impute)  
						writeLines(as.character(GD[1:(ns-1)]), fileNumCon, sep=sep)
						writeLines(as.character(GD[ns]), fileNumCon, sep="\n")
					}#enf of inFIle
				} #end iteration within file
		  
				#Close HMP file
				close.connection(fileHMPCon)
			} #end iteration among files

			#Close GD and GM file
			close.connection(fileNumCon)
			close.connection(fileMapCon)
		}
		if(priority == "speed"){
			fileNumCon<-file(description=paste(out, ".Numeric.txt", sep=""), open="w")
			fileMapCon<-file(description=paste(out, ".map", sep=""), open="w")
			#GM header
			# writeLines("SNP", fileMapCon, sep=sep)
			# writeLines("Chrom", fileMapCon, sep=sep)
			# writeLines("BP", fileMapCon, sep="\n")
			close.connection(fileNumCon)
			close.connection(fileMapCon)
			for(theFile in 1: nFile){
				myFile <- read.table(hfile[theFile], stringsAsFactors=FALSE, colClasses="character", sep=sep, head=FALSE, skip=1)
				nM <- nrow(myFile)
				write.table(myFile[, c(1, 3, 4)], paste(out, ".map", sep=""), append=TRUE, col.names=FALSE, row.names=FALSE, quote=FALSE, sep=sep)
				myFile <- myFile[, -c(1:11)];gc()
				myGDx <- apply(myFile, 1, function(x) KAML.Num(x, impute=SNP.impute))
				myGDx <- t(myGDx)
				write.table(myGDx, paste(out, ".Numeric.txt", sep=""), append=TRUE, col.names=FALSE, row.names=FALSE, quote=FALSE, sep=sep)
				rm(myFile);rm(myGDx);gc()
				cat(paste(" File: ", hfile[theFile], "; Total markers:", nM, " finished!\n", sep=""))
			}
		}
		cat(" Preparation for NUMERIC data is done!\n")
	}

	#Transfer genotype data to .desc, .bin files
	if((!is.null(numfile))|(!is.null(hfile))|(!is.null(vfile))){
		if(is.null(numfile)){
			numfile <- paste(out, ".Numeric.txt", sep="")
			MAP <- read.table(paste(out, ".map", sep=""), stringsAsFactors=FALSE, head=FALSE, sep=sep)
		}else{
			MAP <- read.table(mapfile, head=FALSE, sep=sep, stringsAsFactors=FALSE)
			#write.table(MAP, paste(out, ".map", sep=""), row.names=FALSE, col.names=FALSE, sep=sep, quote=FALSE)
		}
		nmarkers <- nrow(MAP); rm(MAP); gc()
		fileGenoCon <- file(description=numfile, open="r")
		tt2 <-readLines(fileGenoCon, n=1)
		tt2 <- unlist(strsplit(x=tt2, split=sep, fixed=TRUE))
		ns <- length(tt2)
		cat(paste(" Number of Ind: ", ns, "; Number of markers: ", nmarkers, ";\n", sep=""))
		cat(" Output BIG genotype...\n")
		close.connection(fileGenoCon)
		bck <- paste(out, ".geno.bin", sep="")
		dsc <- paste(out, ".geno.desc", sep="")
		if(priority == "memory"){
			myGeno.backed<-big.matrix(nmarkers, ns, type="char",
				backingfile=bck, descriptorfile=dsc)
			#Initialization for iteration within file
			inGENOFile=TRUE
			i=0
			#printN <- unlist(strsplit(x=as.character(nmarkers), split="", fixed=TRUE))
			#printIndex <- seq(0, (as.numeric(printN[1]) + 1) * 10^(length(printN)), 1000)[-1]
			
			#Iteration within file		
			fileGenoCon <- file(description=numfile, open="r")
			while(inGENOFile){
				i=i + maxLine
				tt<-readLines(fileGenoCon, n=maxLine)
				if(i >= nmarkers){
					i <- nmarkers
				}
				# if(sum(i >= printIndex )!=0){
					# printIndex <- printIndex[printIndex > i]
					# print(paste("Number of Markers Written into BIG File: ", i, sep=""))
				# }
				tt<-do.call(rbind, strsplit(x=tt, split=sep, fixed=TRUE))
				if(!is.null(tt) && sum(is.na(tt)) != 0)	stop("'NA' is not allowed in Numeric genotype!")
				nn <- nrow(tt)
				#Identify end of file
				if(is.null(tt[1])) inGENOFile=FALSE
				if(inGENOFile){
					KAML.Bar(i=i, n=nmarkers, fixed.points=FALSE)
					if(i == nmarkers)	cat("\n")
					if(i == nmarkers){
						myGeno.backed [(i-nn + 1):i, ] = tt; rm(tt); gc()
					}else{
						myGeno.backed [(i-maxLine + 1):i, ] = tt; rm(tt); gc()
					}
				}else{
					geno.flush <- flush(myGeno.backed)
				}
			}
			if(!geno.flush){
				stop("flush failed")
			}
			KAML.Impute(myGeno.backed)
			rm(myGeno.backed)
			
			#Close GENO file
			close.connection(fileGenoCon)
		}
		if(priority == "speed"){
			myGeno <- read.big.matrix(numfile, type="char", head=FALSE, sep=sep,
				backingfile=bck, descriptorfile=dsc)
			KAML.Impute(myGeno)
			rm("myGeno")
		}
		cat(" Preparation for BIG data is done!\n")
	}
	gc()
	cat(" KAML data prepration accomplished successfully!\n")
}

KAML.Num <-
function(
	x, impute="Middle"
)
{
#--------------------------------------------------------------------------------------------------------#
# Object: transform ATCG into numeric genotype															 #
# 	 																									 #
# Input:	 																							 #
# x: a vector contains "ATCG"	 																		 #
# impute: "Left", "Middle", "Right"								 										 #
#--------------------------------------------------------------------------------------------------------#
	#replace missing allels
	x[x=="XX"]="N"
	x[x=="--"]="N"
	x[x=="++"]="N"
	x[x=="//"]="N"
	x[x=="NN"]="N"
	x[x=="00"]="N"
	
	#replace false allels
    x[x=="CA"]="AC"
	x[x=="GA"]="AG"
	x[x=="TA"]="AT"
	x[x=="GC"]="CG"
	x[x=="TC"]="CT"
	x[x=="TG"]="GT"
	
    n=length(x)
    lev=levels(as.factor(x))
    lev=setdiff(lev,"N")
    len=length(lev)

    #Genotype counts
    count=1:len
    for(i in 1:len){
        count[i]=length(x[(x==lev[i])])
    }
    position=order(count)
    if(len<=1 | len> 3)x=rep(0, length(x))
    if(len==2)x=ifelse(x=="N",NA,ifelse(x==lev[1],0,2))
	if(len==3)x=ifelse(x=="N",NA,ifelse(x==lev[1],0,ifelse(x==lev[3],2,1)))

    #missing data imputation
    if(impute=="Middle") {x[is.na(x)]=1 }
    if(len==3){
        if(impute=="Minor")  {x[is.na(x)]=position[1]-1}
        if(impute=="Major")  {x[is.na(x)]=position[len]-1}
    }else{
        if(impute=="Minor")  {x[is.na(x)]=2*(position[1]-1)}
        if(impute=="Major")  {x[is.na(x)]=2*(position[len]-1)}
    }
    return(x)
}

KAML.QTN.sel <-
function(
	COR, QTN.list
)
{
#--------------------------------------------------------------------------------------------------------#
# Object: To select the QTNs that increase the predicting accuracy										 #
# 	 																									 #
# Input:	 																							 #
# COR: the vector of accuracy	 																		 #
# QTN.list: a list for storing the QTN							 										 #
#--------------------------------------------------------------------------------------------------------#
	max.pos <-which.max(COR)
	if(max.pos == 1){
		qtn.inc=NULL;qtn.dec=NULL
	}else
	{
		cor.new <- COR[1:max.pos]
		qtn.new <- QTN.list[1:(max.pos-1)]
		cor.a <- rep(1, length(cor.new))
		cor.a[1:(length(cor.a)-1)] <- cor.new[2:length(cor.new)]
		
		cor.sum <- sum((cor.a-cor.new)<0)
		if(cor.sum == 0){
			if(!is.list(qtn.new)){
				qtn.inc=qtn.new[1:(max.pos-1)]
			}else{
				qtn.inc=qtn.new[[max.pos-1]]
			}
			qtn.dec=NULL
		}else
		{
			qtn.dec=NULL
			
			#find the QTN that increase accuracy over the previous QTN and GBLUP
			index <- ((cor.a-cor.new)<0) + (cor.a<cor.new[1])
			#index <- ((cor.a-cor.new)<=0)
			dec.index <- which(index != 0)
			if(is.list(qtn.new)){
				for(d in dec.index){
					if(d == 1){
						qtn.dec=c(qtn.dec, qtn.new[[d]])
					}else
					{
						qtn.dec=c(qtn.dec, setdiff(qtn.new[[d]], qtn.new[[d-1]]))
					}
				}
				qtn.inc=setdiff(qtn.new[[max.pos-1]], qtn.dec)
			}else{
				qtn.dec <- qtn.new[dec.index]
				qtn.inc <- setdiff(qtn.new, qtn.dec)
			}
		}
	}
	return(list(qtn.inc=qtn.inc, qtn.dec=qtn.dec))
}

KAML.CrossV <- 
function(
	phe, geno, K=NULL, CV=NULL, GWAS.model=NULL, qtn.model="SR", BF.threshold=NULL, vc.method="emma", Top.num=NULL, Top.perc=NULL, max.nQTN=TRUE, SNP.filter=0.5,
	cor.threshold=0.99, count.threshold=0.9, sample.num=1, crv.num=5, cpu=1, theSeed=NULL, prior.QTN=NULL, K.method="center", binary=FALSE, priority="speed",
	bin.size=1000000, amplify=NULL, bisection.loop=5, ref.gwas=FALSE, prior.model="QTN+K", GWAS.npc=3
)
{
#--------------------------------------------------------------------------------------------------------#
# Object: To perform cross-validation in reference population											 #
# 	 																									 #
# Input:	 																							 #
# phe: Phenotype, a value vector(NA is allowed, only non-NA individuals will be used for reference)		 #
# geno: Genotype in numeric format, pure 0, 1, 2 matrix; m * n, m is marker size, n is population size	 #
# (Note: both matrix or big.matrix are acceptable)														 #
# K: Kinship for all individuals(row and column orders must be the same with phenotype)					 #
# CV: covariance, design matrix(n * x) for the fixed effects(CV must contain a column of 1's)			 #
# map: the map for genotype; three columns: snp, chrom, position										 #
# GWAS.model: which model will be used for GWAS(only "GLM" and "MLM" can be selected presently)			 #
# vc.method: method for variance components estimation("emma" or "gemma")								 #
# Top.num: a number, a subset of top SNPs for each iteration are used as covariants						 #
# Top.perc: a vector, a subset of top SNPs for each iteration are amplified when calculating KINSHIP	 #
# max.nQTN: whether limits the max number of Top.index													 #
# SNP.filter: the P-values of T.test for SNPs which below this threshold will be deleted				 #
# cor.threshold: if the top SNP which is used as covariant is in high correlation with others,			 #
# 			it will be deleted																			 #
# sample.num: the sample number of cross-validation														 #
# crv.num: the cross number of cross-validation															 #
# count.threshold: if the count of selected SNP for all iteration >= sample.num*crv.num*count.threshold, #
# 			than it will be treated as covariance in final predicting model								 # 
# cpu: the number of CPU for calculation																 #
# theSeed: the random seed 																				 #
# prior.QTN: the prior QTNs which will be added as covariants, if provided prior QTNs, 	KAML will not	 #
# 			optimize QTNs and model during cross-validation												 #
# prior.model: the prior MODEL(only "K", "QTN+K", "QTN" can be selected presently)						 #
# K.method: which algorithm will be applied("center", "scale", "vanraden")								 #
# bin.size: the size of each bin																		 #
# amplify: a vector, the base for LOG																	 #
# bisection: whether using bisection algorithm to optimize KINSHIP										 #
# bisection.loop: the max loop(iteration) number of bisection algorithm									 #
# ref.gwas: whether to do GWAS for reference population(if not, KAML will merge all GWAS results of	 #
# 			cross-validation by mean)																	 #
# GWAS.npc=3: the number of PC that will be added as covariance to control population structure			 #
#--------------------------------------------------------------------------------------------------------#
	#make sure the type of system(windows, mac, linux)
	R.ver <- Sys.info()[['sysname']]
	wind <- R.ver == 'Windows'
	linux <- R.ver == 'Linux'
	mac <- (!linux) & (!wind)
	
	#make sure the type of R(base R, Open R), "checkpoint" package is installed in Open R by default
    #r.open <- "checkpoint" %in% rownames(installed.packages())
    r.open <- !inherits(try(Revo.version, silent=TRUE),"try-error")
	mkl <- Math_cpu_check()
	if(r.open &  mac & cpu > 1)	Sys.setenv("VECLIB_MAXIMUM_THREADS" = "1")
	if(wind)	cpu <- 1

	if(is.null(Top.num) & (is.null(Top.perc))) stop("One of 'Top.num' or 'Top.perc' must be setted")
	if(!qtn.model %in% c("MR", "SR", "BF"))	stop("Please select the right 'qtn.model': 'MR', 'SR', 'BF'")
	
	N.Ind <- length(phe)
	NA.index <- is.na(phe)
	NA.Ind <- sum(NA.index)
	
	n.ref <- N.Ind-NA.Ind
	n.inf <- NA.Ind
	n.marker <- nrow(geno)

	if(!is.null(theSeed)) {set.seed(theSeed)}
	if(!is.null(prior.QTN))	Top.num <- NULL
	if(qtn.model == "BF" & is.null(BF.threshold))	BF.threshold <- 0.05 / n.marker
	
	inf.index <- list()
	P.value <- list()
	if(!is.null(Top.num) | (!is.null(Top.perc) && (length(Top.perc) > 1 | length(amplify) > 1))){
		
		cat(paste(" Cross-validation(", sample.num, "*", crv.num, ") Started...\n", sep=""))
		#get the inference population index for cross-validation
		for(s in 1:sample.num){
			CrossVindex <- sample(which(!NA.index), n.ref)
			for(cv in 1:crv.num){
				if(cv < crv.num){
					cv.index <- ((cv-1) * (n.ref%/%crv.num) + 1):(cv * (n.ref%/%crv.num))
				}else
				{
					cv.index <- ((cv-1) * (n.ref%/%crv.num) + 1):(n.ref)
				}
				inf.index[[(s-1)*crv.num+cv]] <- CrossVindex[cv.index]
			}
		}

		#change GWAS models at here
		mult.run.gwas <- function(i){
			if(r.open)	try(setMKLthreads(mkl), silent=TRUE)
			ref.logic <- (i == sample.num * crv.num + 1)
			# cat(" GWAS of validations is ongoing...(", i,"/", gwas.num, ")\r", sep="")
			if(!ref.logic){
				P.ref <- phe
				P.ref[inf.index[[i]]] <- NA
			}else{
				P.ref <- phe
			}
			if(!ref.logic){
				pri <- paste(" GWAS of validations NO.", i, sep="")
			}else{
				pri <- " GWAS of Total References finished"
			}
			if(GWAS.model == "GLM"){
				# if(!is.null(SNP.filter)){
				# 	GLM.gwas <- KAML.GWAS(phe=P.ref, geno=deepcopy(geno[SNP.index, , drop=FALSE], CV=CV, method="GLM", cpu=cpus, priority=priority, NPC=GWAS.npc, bar.head=pri, bar.tail="", bar.len=0)
				# 	P.value <- vector("numeric", length=n.marker)
				# 	P.value[SNP.index] <- GLM.gwas[, 2]
				# 	P.value[-SNP.index] <- SNP.filter
				# }else{
					GLM.gwas <- KAML.GWAS(phe=P.ref, geno=geno, CV=CV, method="GLM", cpu=cpus, priority=priority, NPC=GWAS.npc, bar.head=pri, bar.tail="", bar.len=1)
					P.value <- GLM.gwas[, 3]
				# }
				rm("GLM.gwas")
			}
			if(GWAS.model == "RR"){
					rr.gwas <- KAML.Mix(phe=P.ref, geno=geno, CV=CV, K=K, vc.method=vc.method)$S
					P.value <- as.vector(geno[, !is.na(P.ref), drop=FALSE] %*% rr.gwas)/nrow(geno)
					P.value <- -(abs(P.value)/mean(abs(P.value)))
					P.value <- 10^(P.value)
					rm("rr.gwas")
			}
			if(GWAS.model == "MLM"){
				# if(!is.null(SNP.filter)){
				# 	MLM.gwas <- KAML.GWAS(phe=P.ref, geno=geno[SNP.index, , drop=FALSE], CV=CV, K=K, method="MLM", vc.method=vc.method, cpu=cpus, NPC=GWAS.npc, bar.head=pri, bar.tail="", bar.len=0)
				# 	P.value <- vector("numeric", length=n.marker)
				# 	P.value[SNP.index] <- MLM.gwas[, 2]
				# 	P.value[-SNP.index] <- SNP.filter
				# }else{
					MLM.gwas <- KAML.GWAS(phe=P.ref, geno=geno, CV=CV, K=K, method="MLM", vc.method=vc.method, cpu=cpus, priority=priority, NPC=GWAS.npc, bar.head=pri, bar.tail="", bar.len=1)
					P.value <- MLM.gwas[, 3]
					rm("MLM.gwas")
				# }
			}
			# if(i == gwas.num)	cat("\n")
			rm(list=c("P.ref")); gc()
			P.value[is.na(P.value)] <- 1
			P.value[P.value == 0] <- min(P.value[P.value != 0])
			return(P.value)
		}
		
		GWAS.model.txt <- ifelse(GWAS.model=="GLM", "GLM(Generalized Linear Model)", ifelse(GWAS.model=="MLM", "MLM(Mixed Linear Model)", "RR(Ridge Regression)"))
		cat(paste(" GWAS with model: ", GWAS.model.txt, "\n", sep=""))
		
		#the SNPs which are lower associated than 'SNP.filter' will be deleted to save time
		P.value.ref <- NULL
		# if(!is.null(SNP.filter)){
		# 	cat(paste(" Filtering SNPs with the threshold: ", SNP.filter, "\n", sep=""))
		# 	GLM.gwas <- KAML.GWAS(phe=phe, geno=geno, CV=CV, method="GLM", NPC=GWAS.npc, cpu=cpu, bar.head=" Filtering finished ", bar.tail="", bar.len=0)
		# 	P.value.ref <- GLM.gwas[, 2]
		# 	P.value.ref[is.na(P.value.ref)] <- 1
		# 	SNP.index <- which(P.value.ref < SNP.filter)
		# 	cat(paste(" Number of SNPs deleted: ", n.marker-length(SNP.index), "\n", sep=""))
		# 	if(GWAS.model == "GLM" & ref.gwas){
		# 		rm("GLM.gwas"); gc()
		# 	}else{
		# 		P.value.ref=NULL; rm("GLM.gwas"); gc()
		# 	}
		# }
		
		if(ref.gwas & is.null(P.value.ref)){
			gwas.num <- sample.num * crv.num + 1
		}else{
			gwas.num <- sample.num * crv.num
		}

		# if(wind & cpu > 1){
		# 	cpus <- 1
		# 	cat(" Multi-process of GWAS started...\n")
		# 	max.cpu <- min(cpu, gwas.num)
		# 	cl <- makeCluster(max.cpu, outfile = "Loop.log")
		# 	registerDoParallel(cl)
		# 	clusterExport(cl, varlist=c("KAML.GWAS", "KAML.EMMA.REML", "KAML.Mix"))
		# 	P.values <- foreach(x=1:gwas.num, .packages=c("bigmemory", "rfunctions")) %dopar% mult.run.gwas(x)
		# 	stopCluster(cl)
		# 	cat(" Multi-process done!\n")
		# }else{
			cpus <- cpu
			P.values <- lapply(1:gwas.num, mult.run.gwas)
		# }
		
		#----------debug-----------------#
		#wrt.p <- do.call(cbind, P.values)
		#write.csv(wrt.p,"cv.p.csv")
		#--------------------------------#

		P.value <- P.values[1:(sample.num * crv.num)]
        if(ref.gwas && is.null(P.value.ref)){
            P.value.ref <- unlist(P.values[sample.num * crv.num + 1])
        }else if(ref.gwas && !is.null(P.value.ref)){
			P.value.ref <- P.value.ref
		}else{
			cat(" Merge the GWAS of Cross-validations by: mean\n")
			P.value.ref <- apply(matrix(unlist(P.value), n.marker), 1, mean)
        }
		rm("P.values"); gc()

		#----------------------------------debug-----------------------------------------#
		# library(CMplot)
		 # GWAS.res <- cbind(map, P.value.ref)
		 # GWAS.res[,1] <- 1:nrow(GWAS.res)
		 # CMplot(GWAS.res, plot.type="m",threshold=0.05,file="pdf",col=c("black","orange"))
		 # stop()
		#--------------------------------------------------------------------------------#
			
		if(binary){
			mytext <- "AUC"
			mytext1 <- "AUC(Area Under the Curve)"
		}else{
			mytext <- "COR"
			mytext1 <- "COR(Correlation)"
		}
		cat(paste(" Benchmark of Performance:", mytext1, "\r\n"))
	}else{
		cat(paste(" Cross-validation(References only) Started...\n", sep=""))
		cat(paste(" GWAS with model: ", GWAS.model, "...\n", sep=""))
		P.ref <- phe
		if(GWAS.model == "GLM"){
			GLM.gwas <- KAML.GWAS(phe=P.ref, geno=geno, CV=CV, method="GLM", NPC=GWAS.npc, cpu=cpu, priority=priority, bar.head=" GWAS of References", bar.tail="", bar.len=1)
			P.value.ref <- GLM.gwas[, 3]
			rm("GLM.gwas")
		}
		if(GWAS.model == "RR"){
			rr.gwas <- KAML.Mix(phe=P.ref, geno=geno, CV=CV, K=K, vc.method=vc.method)$S
			P.value <- as.vector(geno[, !is.na(P.ref), drop=FALSE] %*% rr.gwas)/nrow(geno)
			P.value <- -(abs(P.value)/mean(abs(P.value)))
			P.value <- 10^(P.value)
			rm("rr.gwas")
		}
		if(GWAS.model == "MLM"){
			MLM.gwas <- KAML.GWAS(phe=P.ref, geno=geno, CV=CV, K=K, method="MLM", vc.method=vc.method, NPC=GWAS.npc, cpu=cpu, priority=priority, bar.head=" GWAS of References", bar.tail="", bar.len=1)
			P.value.ref <- MLM.gwas[, 3]
			rm("MLM.gwas")
		}
		P.value.ref[is.na(P.value.ref)] <- 1
		rm(list=c("P.ref")); gc()
	}

	sam <- sample.num
	#sample.num <- 1
	#cpu <- 1
	
	#Cross-validation to optimize pseudo QTN
	if(is.null(Top.num)){
		if(!is.null(prior.QTN)){
			cross.QTN = prior.QTN
			cross.model=prior.model
			cat(" The provided MODEL:\n")
			cat(" ", cross.model, "\n")
			cat(" The provided QTNs:\n")
			cat(" ", cross.QTN, "\n")
		}else{
			cat(" No QTNs and MODEL optimization\n")
			cross.QTN = NULL
			cross.model="K"
		}
	}else
	{
		pseudo.QTN.k <-list()
		pseudo.QTN.qtn <-list()
		bin.qtn <- list()
		bin.sel <- list()
		inc.QTN.store.k <- NULL
		inc.QTN.store.qtn <- NULL
		cor.pseudo.k <- NULL
		cor.pseudo.qtn <- NULL
		
		Top.index <- c(0, 1 : Top.num)
		uni <- (n.ref) %/% crv.num
		
		#set the upper bound of pseudo QTNs
		if(max.nQTN == TRUE){
			max.nQTN.v <- round(sqrt(n.ref) / sqrt(log10(n.ref)))
			Top.index[Top.index >= max.nQTN.v] <- max.nQTN.v
			Top.index <- sort(unique(Top.index))
		}
		Top.index <- Top.index[Top.index < (crv.num-1) * uni]
		
		#----------debug----------#
		#print(Top.index)
		#-------------------------#
		
		cat(" Pick up pseudo QTNs at LD threshold(", cor.threshold, ") as following: \n", sep="")
		
		for(g in 1: length(P.value[1:(sample.num * crv.num)])){
			P.value.index <- order(P.value[[g]])
			bin.qtn[[g]] <- P.value.index[1]
			reps <- 1
			repeat({
				reps <- reps + 1
				p_cor <- apply(geno[bin.qtn[[g]], , drop=FALSE], 1, function(x){abs(cor(x, geno[P.value.index[reps], ]))})
				if(sum(p_cor > cor.threshold) == 0)	bin.qtn[[g]] <- c(bin.qtn[[g]], P.value.index[reps])
				if(length(bin.qtn[[g]]) == max(Top.index) | reps == length(P.value.index))	break()
			})

			#----------------debug----------------#
			 cat(" ", as.vector(bin.qtn[[g]]), "\n")
			 #cat(" ", bin.sel[[g]], "\n")
			#-------------------------------------#	
		}
		rm(list=c("P.value.index"));gc()
		cat(" QTNs are confirmed!\n")
		cat(" Optimizing pseudo QTNs and MODEL...\n")
		bin.qtn.unique <- unique(unlist(bin.qtn))
		bin.qtn.unique <- bin.qtn.unique[order(P.value.ref[bin.qtn.unique])]
		bin.count <- table(unlist(bin.qtn))
		bin.names <- names(bin.count)
		bin.count <- bin.count[match(bin.qtn.unique, bin.names)]
		
		qtn.cor.matrix <- abs(cor(t(geno[bin.qtn.unique, , drop=FALSE])))
		diag(qtn.cor.matrix) <- 0
		replace.index <- which(qtn.cor.matrix > cor.threshold, arr.ind=TRUE)
		if(nrow(replace.index) != 0){
			replace.index <- replace.index[replace.index[, 2] > replace.index[, 1], , drop=FALSE]
			replace.index <- replace.index[order(replace.index[, 1]), , drop=FALSE]
			for(bin.q in 1:nrow(replace.index)){
				bin.qtn <- lapply(bin.qtn, function(x){
					x[x==bin.qtn.unique[replace.index[bin.q, 2]]] <- bin.qtn.unique[replace.index[bin.q, 1]]; return(x)}
				)
			}
		}
		
		bin.count <- table(unlist(bin.qtn))
		bin.names <- names(bin.count)
		bin.pseudo <- as.numeric(bin.names[bin.count >= floor(sample.num * crv.num * count.threshold)])
		if(length(bin.pseudo) == 0){
			cat(" Prior pseudo QTNs:\n")
			cat(" NULL\n")
			pseudo.QTNs <- NULL
		}else{
			pseudo.QTNs <- bin.pseudo[order(P.value.ref[bin.pseudo])]
			cat(" Prior pseudo QTNs:\n")
			cat(" ", pseudo.QTNs, "\n")
		}
		
		Top.index <- Top.index[Top.index <= length(pseudo.QTNs)]

		if(is.null(pseudo.QTNs)){
			cross.QTN <- NULL
			cross.model <- "K"
		}else
		{
			qtn.model.txt <- ifelse(qtn.model == "MR", "Multiple Regression", ifelse(qtn.model == "SR", "Step Regression", "Top picked"))
			cat(" Choosed algorithm:", qtn.model.txt, "\n")
			
			#do parallel
			mr.mult.run <- function(j, math.cpu=NULL){
			
				#when use microsoft r open, users can set the number of cpu for math calculation at each process
				#if(cpu>1 & (wind | r.open)) try(setMKLthreads(math.cpu), silent=TRUE)
				
				top <- j %% length(Top.index)
				if(top == 0){
					top <- length(Top.index)
					row.index <- j%/% length(Top.index)
				}else
				{
					row.index <- j%/% length(Top.index) + 1
				}
				mult.inf.index <- inf.index[[row.index]]
				mult.ref.index <- !NA.index
				mult.ref.index[mult.inf.index] <- FALSE

				P.ref <- phe
				P.ref[!mult.ref.index] <- NA
				P.inf <- phe[mult.inf.index]

				if((Top.index[top]) == 0){
					gblup <- KAML.Mix(phe=P.ref, K=K, CV=CV, vc.method=vc.method)
					gebv <- (CV %*% as.matrix(gblup$beta) + gblup$ebv)[mult.inf.index]
					if(binary){
						acc.k <- KAML.stas.cal(P.inf, gebv, type="auc")
						acc.qtn <- KAML.stas.cal(P.inf, gebv, type="auc")
					}else{
						acc.k <- KAML.stas.cal(P.inf, gebv, type="cor")
						acc.qtn <- KAML.stas.cal(P.inf, gebv, type="cor")
					}
					QTN.store <- NULL
					#cat(paste(" Cross-validation NO.", row.index, "; NQTN=", Top.index[top], "; Model='K'", "; ",mytext ,"=", round(acc.k, 4), paste(rep(" ", 5), collapse=""), "\r", sep=""))
				}else
				{
					myqtn <- pseudo.QTNs[1:(Top.index[top])]
					QTN.cv <- t(geno[myqtn, mult.ref.index, drop=FALSE])

					#calculate GLM model with QTNs
					qtn.eff <- KAML.GLM(y=P.ref[mult.ref.index], X=CV[mult.ref.index, ], qtn.matrix=QTN.cv)$QTN.eff
					if(length(myqtn) == 1){
						gebv <- cbind(CV[mult.inf.index, ], geno[myqtn, mult.inf.index]) %*% as.matrix(qtn.eff)
						if(binary){
							acc.qtn <- KAML.stas.cal(P.inf, gebv, type="auc")
						}else{
							acc.qtn <- KAML.stas.cal(P.inf, gebv, type="cor")
						}
					}else{
						gebv <- cbind(CV[mult.inf.index, ], t(geno[myqtn, mult.inf.index])) %*% as.matrix(qtn.eff)
						if(binary){
							acc.qtn <- KAML.stas.cal(P.inf, gebv, type="auc")
						}else{
							acc.qtn <- KAML.stas.cal(P.inf, gebv, type="cor")
						}
					}
					
					#cat(paste(" Cross-validation NO.", row.index, "; NQTN=", Top.index[top], "; Model='QTN'", "; ", mytext, "=", round(acc.qtn, 4), paste(rep(" ", 5), collapse=""), "\r", sep=""))
					
					#calculate MLM model with QTNs
					KAML.cv <- t(geno[myqtn, , drop=FALSE])
					KAML.cv <- cbind(CV, KAML.cv)
					gblup <- KAML.Mix(phe=P.ref, vc.method=vc.method, CV=KAML.cv, K=K)
					gebv <- (gblup$ebv + KAML.cv %*% as.matrix(gblup$beta))[mult.inf.index]
					if(binary){
						acc.k <- KAML.stas.cal(P.inf, gebv, type="auc")
					}else{
						acc.k <- KAML.stas.cal(P.inf, gebv, type="cor")
					}
					#textx <- paste(" Cross-validation NO.", row.index, "; NQTN=", Top.index[top], "; Model='QTN+K'", "; ", mytext, "=", round(acc.k, 4), sep="")
					#cat( textx, paste(rep(" ", 5), collapse=""), "\r")
				}
				print.f(j)
				rm(list=c("P.ref", "P.inf", "gblup", "gebv")); gc()
				acc.qtn[is.na(acc.qtn)] <- 0
				acc.k[is.na(acc.k)] <- 0
				return(list(acc.k=acc.k, acc.qtn=acc.qtn))
			}

			sr.mult.run <- function(j, match.cpu=NULL){				
				mult.inf.index <- inf.index[[j]]
				mult.ref.index <- !NA.index
				mult.ref.index[mult.inf.index] <- FALSE

				P.ref <- phe
				P.ref[!mult.ref.index] <- NA
				P.inf <- phe[mult.inf.index]
				if(is.null(myqtn)){
					gblup <- KAML.Mix(phe=P.ref, K=K, CV=CV, vc.method=vc.method)
					gebv <- (CV %*% as.matrix(gblup$beta) + gblup$ebv)[mult.inf.index]
					if(binary){
						acc.k <- KAML.stas.cal(P.inf, gebv, type="auc")
						acc.qtn <- KAML.stas.cal(P.inf, gebv, type="auc")
					}else{
						acc.k <- KAML.stas.cal(P.inf, gebv, type="cor")
						acc.qtn <- KAML.stas.cal(P.inf, gebv, type="cor")
					}
					#cat(paste(" Cross-validation NO.", row.index, "; NQTN=", Top.index[top], "; Model='K'", "; ",mytext ,"=", round(acc.k, 4), paste(rep(" ", 5), collapse=""), "\r", sep=""))
				}else
				{
					glm.myqtn <- myqtn[glm.qtn.index[1 : Top.index[top]]]
					QTN.cv <- t(geno[glm.myqtn, mult.ref.index, drop=FALSE])

					#calculate GLM model with QTNs
					qtn.eff <- KAML.GLM(y=P.ref[mult.ref.index], X=CV[mult.ref.index, , drop=FALSE], qtn.matrix=QTN.cv)$QTN.eff
					gebv <- cbind(CV[mult.inf.index, , drop=FALSE], t(geno[glm.myqtn, mult.inf.index, drop=FALSE])) %*% as.matrix(qtn.eff)
					if(binary){
						acc.qtn <- KAML.stas.cal(P.inf, gebv, type="auc")
					}else{
						acc.qtn <- KAML.stas.cal(P.inf, gebv, type="cor")
					}

					#cat(paste(" Cross-validation NO.", row.index, "; NQTN=", Top.index[top], "; Model='QTN'", "; ", mytext, "=", round(acc.qtn, 4), paste(rep(" ", 5), collapse=""), "\r", sep=""))
					
					#calculate MLM model with QTNs
					mlm.myqtn <- myqtn[mlm.qtn.index[1 : Top.index[top]]]					
					KAML.cv <- t(geno[mlm.myqtn, , drop=FALSE])
					KAML.cv <- cbind(CV, KAML.cv)
					gblup <- KAML.Mix(phe=P.ref, vc.method=vc.method, CV=KAML.cv, K=K)
					gebv <- (gblup$ebv + KAML.cv %*% as.matrix(gblup$beta))[mult.inf.index]
					if(binary){
						acc.k <- KAML.stas.cal(P.inf, gebv, type="auc")
					}else{
						acc.k <- KAML.stas.cal(P.inf, gebv, type="cor")
					}
				}
				
				print.f(j + (top - 1) * (sample.num * crv.num))

				return(list(acc.k=acc.k, acc.qtn=acc.qtn))
			}
			
			iterationN <- ((sample.num * crv.num) * length(Top.index))
			if(qtn.model == "MR"){mult.run <- mr.mult.run; loop <- ((sample.num * crv.num) * length(Top.index))}
			if(qtn.model == "SR"){mult.run <- sr.mult.run; loop <- (sample.num * crv.num)}
			
			if(qtn.model != "BF"){
				cat(paste(" Total iteration number:", iterationN))
				if(cpu == 1){
					print.f <- function(i){KAML.Bar(i=i, n=iterationN, type="type1", fixed.points=FALSE, symbol.len=0, symbol.head=" Cross-validation Finished_", symbol.tail="")}
					if(qtn.model == "MR"){
						mult.res <- lapply(1:loop, mult.run); gc()
						res <- do.call(rbind, mult.res)
						cor.store.k <- matrix(unlist(res[, 1]), length(Top.index))		
						cor.store.qtn <- matrix(unlist(res[, 2]), length(Top.index))
						
						#----------debug----------#
						# print(cor.store.k)
						# print(cor.store.qtn)
						#-------------------------#
						
						cor.store.k_mean <- rowMeans(cor.store.k)
						cor.store.qtn_mean <- rowMeans(cor.store.qtn)
						
						#----------debug----------#
						# print(pseudo.QTNs)
						# print(cor.store.k_mean)
						#-------------------------#					
						
						#Select the qtns which can increase the accuracy for each validation
						inc.QTN.store.k <- KAML.QTN.sel(COR=cor.store.k_mean, QTN.list=pseudo.QTNs)$qtn.inc
						inc.QTN.store.qtn <- KAML.QTN.sel(COR=cor.store.qtn_mean, QTN.list=pseudo.QTNs)$qtn.inc
						
						#-----------debug-----------#
						# print(inc.QTN.store.k)
						#---------------------------#
						
						k.inc.index <- (cor.store.k[2:nrow(cor.store.k), , drop=FALSE] - cor.store.k[1:(nrow(cor.store.k)-1), , drop=FALSE]) > 0
						k.inc.index <- apply(k.inc.index, 1, sum)
						inc.QTN.store.k <- intersect(inc.QTN.store.k, pseudo.QTNs[k.inc.index >= floor(sample.num * crv.num * count.threshold)])
						qtn.inc.index <- (cor.store.qtn[2:nrow(cor.store.qtn), , drop=FALSE] - cor.store.qtn[1:(nrow(cor.store.qtn)-1), , drop=FALSE]) > 0
						qtn.inc.index <- apply(qtn.inc.index, 1, sum)
						inc.QTN.store.qtn <- intersect(inc.QTN.store.qtn, pseudo.QTNs[qtn.inc.index >= floor(sample.num * crv.num * count.threshold)])

					}
					if(qtn.model == "SR"){
						glm.cor.store <- NULL
						mlm.cor.store <- NULL
						glm.logic.index <- list()
						mlm.logic.index <- list()
						glm.qtn.index <- rep(TRUE, max(Top.index))
						mlm.qtn.index <- rep(TRUE, max(Top.index))
						for(top in 1:length(Top.index)){
							if(top == 1){
								myqtn <- NULL
							}else{
								myqtn <- pseudo.QTNs[1:(Top.index[top])]
							}
							mult.res <- lapply(1:loop, mult.run, mc.cores = cpus)
							res <- do.call(rbind, mult.res)
							cor.store.k <- matrix(unlist(res[, 1]), length(Top.index))
							mlm.logic.index[[top]] <- cor.store.k
							cor.store.qtn <- matrix(unlist(res[, 2]), length(Top.index))
							glm.logic.index[[top]] <- cor.store.qtn
							glm.cor.store[top] <- mean(cor.store.qtn)
							mlm.cor.store[top] <- mean(cor.store.k)
							if(top > 2){
								if((glm.cor.store[top] < glm.cor.store[top - 1]) & sum(glm.logic.index[[top]]-glm.logic.index[[top-1]]) >= floor(sample.num * crv.num * count.threshold)){
									glm.cor.store[top] <- glm.cor.store[top - 1]
									glm.qtn.index[(Top.index[top - 1] + 1) : Top.index[top]] <- FALSE
								}
							}
							if(top > 1){
								if((mlm.cor.store[top] < mlm.cor.store[top - 1]) & sum(mlm.logic.index[[top]]-mlm.logic.index[[top-1]]) >= floor(sample.num * crv.num * count.threshold)){
									mlm.cor.store[top] <- mlm.cor.store[top - 1]
									mlm.qtn.index[(Top.index[top - 1] + 1) : Top.index[top]] <- FALSE
								}
							}
						}
						inc.QTN.store.qtn <- pseudo.QTNs[glm.qtn.index]
						inc.QTN.store.k <- pseudo.QTNs[mlm.qtn.index]
					}
				}else
				{
					#1：mclapply doesn't work at windows system
					#2: there are always some errors when use multi-process in microsoft r open
					
					#get the user-specific cpu number
					if(wind){
						print.f <- function(i){KAML.Bar(i=i, n=iterationN, type="type1", fixed.points=FALSE, symbol.len=0, symbol.head=" Cross-validation Finished_", symbol.tail="")}
						#print(paste("Cross-validation NO.", paste(c(1:(sample.num * crv.num)), collapse=", "), sep=""))
						cat(" Multi-process started...\n")
						cat(" (Note: There needs to wait some time! See iteration details in 'Loop.log')\n")
						cl <- makeCluster(cpu, outfile = "Loop.log")
						registerDoParallel(cl)
						clusterExport(cl, varlist=c("KAML.Mix", "KAML.EMMA.REML", "KAML.HE", "KAML.GEMMA.VC",
							"KAML.stas.cal", "KAML.Kin", "KAML.QTN.rm", "KAML.GLM"))
						mult.res <- foreach(x=1:loop,
						.packages=c("bigmemory", "rfunctions")) %dopar% mult.run(x, math.cpu=mkl)
						if(is.null(Top.perc))	stopCluster(cl)
						cat(" Multi-process done!\n")
					}else{
						# print.f <- function(i){writeBin(1, tmpf)}
						# KAML.Bar(n=iterationN, type="type2", tmp.file=tmpf, fixed.points=FALSE, symbol.len=0, symbol.head=" Cross-validation Finished_", symbol.tail="")

						if(r.open)	setMKLthreads(1)
						cpus <- ifelse(loop > cpu, cpu, loop)
						if(qtn.model == "MR"){
							tmpf.name <- tempfile()
							tmpf <- fifo(tmpf.name, open="w+b", blocking=TRUE)
							writeBin(0, tmpf)
							print.f <- function(i){KAML.Bar(n=iterationN, type="type3", tmp.file=tmpf, fixed.points=FALSE, symbol.len=0, symbol.head=" Cross-validation Finished_", symbol.tail="")}
							mult.res <- parallel::mclapply(1:loop, mult.run, mc.cores = cpus)
							close(tmpf); unlink(tmpf.name)
							res <- do.call(rbind, mult.res)
							cor.store.k <- matrix(unlist(res[, 1]), length(Top.index))		
							cor.store.qtn <- matrix(unlist(res[, 2]), length(Top.index))
							
							#----------debug----------#
							# print(cor.store.k)
							# print(cor.store.qtn)
							#-------------------------#
							
							cor.store.k_mean <- rowMeans(cor.store.k)
							cor.store.qtn_mean <- rowMeans(cor.store.qtn)
							
							#----------debug----------#
							# print(pseudo.QTNs)
							# print(cor.store.k_mean)
							#-------------------------#					
							
							#Select the qtns which can increase the accuracy for each validation
							inc.QTN.store.k <- KAML.QTN.sel(COR=cor.store.k_mean, QTN.list=pseudo.QTNs)$qtn.inc
							inc.QTN.store.qtn <- KAML.QTN.sel(COR=cor.store.qtn_mean, QTN.list=pseudo.QTNs)$qtn.inc
							
							#-----------debug-----------#
							# print(inc.QTN.store.k)
							#---------------------------#
							
							k.inc.index <- (cor.store.k[2:nrow(cor.store.k), , drop=FALSE] - cor.store.k[1:(nrow(cor.store.k)-1), , drop=FALSE]) > 0
							k.inc.index <- apply(k.inc.index, 1, sum)
							inc.QTN.store.k <- intersect(inc.QTN.store.k, pseudo.QTNs[k.inc.index >= floor(sample.num * crv.num * count.threshold)])
							qtn.inc.index <- (cor.store.qtn[2:nrow(cor.store.qtn), , drop=FALSE] - cor.store.qtn[1:(nrow(cor.store.qtn)-1), , drop=FALSE]) > 0
							qtn.inc.index <- apply(qtn.inc.index, 1, sum)
							inc.QTN.store.qtn <- intersect(inc.QTN.store.qtn, pseudo.QTNs[qtn.inc.index >= floor(sample.num * crv.num * count.threshold)])
						
						}
						if(qtn.model == "SR"){
							glm.cor.store <- NULL
							mlm.cor.store <- NULL
							glm.logic.index <- list()
							mlm.logic.index <- list()
							glm.qtn.index <- rep(TRUE, max(Top.index))
							mlm.qtn.index <- rep(TRUE, max(Top.index))
							for(top in 1:length(Top.index)){
								if(top == 1){
									myqtn <- NULL
								}else{
									myqtn <- pseudo.QTNs[1:(Top.index[top])]
								}
								tmpf.name <- tempfile()
								tmpf <- fifo(tmpf.name, open="w+b", blocking=TRUE)
								writeBin(0, tmpf)
								print.f <- function(i){KAML.Bar(n=iterationN, type="type3", tmp.file=tmpf, fixed.points=FALSE, symbol.len=0, symbol.head=" Cross-validation Finished_", symbol.tail="")}
								mult.res <- parallel::mclapply(1:loop, mult.run, mc.cores = cpus)
								close(tmpf); unlink(tmpf.name)
								res <- do.call(rbind, mult.res)
								cor.store.k <- matrix(unlist(res[, 1]), length(Top.index))
								mlm.logic.index[[top]] <- cor.store.k
								cor.store.qtn <- matrix(unlist(res[, 2]), length(Top.index))
								glm.logic.index[[top]] <- cor.store.qtn
								glm.cor.store[top] <- mean(cor.store.qtn)
								mlm.cor.store[top] <- mean(cor.store.k)
								if(top > 2){
									if((glm.cor.store[top] < glm.cor.store[top - 1]) & sum(glm.logic.index[[top]]-glm.logic.index[[top-1]]) >= floor(sample.num * crv.num * count.threshold)){
										glm.cor.store[top] <- glm.cor.store[top - 1]
										glm.qtn.index[(Top.index[top - 1] + 1) : Top.index[top]] <- FALSE
									}
								}
								if(top > 1){
									if((mlm.cor.store[top] < mlm.cor.store[top - 1]) & sum(mlm.logic.index[[top]]-mlm.logic.index[[top-1]]) >= floor(sample.num * crv.num * count.threshold)){
										mlm.cor.store[top] <- mlm.cor.store[top - 1]
										mlm.qtn.index[(Top.index[top - 1] + 1) : Top.index[top]] <- FALSE
									}
								}
							}
							inc.QTN.store.qtn <- pseudo.QTNs[glm.qtn.index]
							inc.QTN.store.k <- pseudo.QTNs[mlm.qtn.index]
							cor.store.qtn <- max(glm.cor.store)
							cor.store.k <- max(mlm.cor.store)
						}
						if(r.open)	setMKLthreads(mkl)
					}
					gc()
				}
	
				if(length(inc.QTN.store.k) == 0)	inc.QTN.store.k <- NULL
				if(length(inc.QTN.store.qtn) == 0)	inc.QTN.store.qtn <- NULL
				
				#-----------debug-----------#
				# print(inc.QTN.store.k)
				#---------------------------#

			}else{
				BinQTN <- pseudo.QTNs
				inc.QTN.store.k <- BinQTN[P.value.ref[BinQTN] < BF.threshold]
				inc.QTN.store.qtn <- NULL
				model.index <- 2
			}
						
			if(qtn.model == "MR"){
				cor.gblup <- cor.store.k[1, ]
				cor.mean.k <- rowMeans(cor.store.k)
				cor.mean.qtn <- rowMeans(cor.store.qtn)
				cor.mean.k.max <- max(cor.mean.k[-1])
				cor.mean.qtn.max <- max(cor.mean.qtn[-1])
				model.index <- which.max(c(cor.mean.k[1], cor.mean.k.max, cor.mean.qtn.max))
				cor.mean <- rbind(cor.mean.qtn, cor.mean.k)
				rownames(cor.mean)=c("QTN", "QTN+K")
				colnames(cor.mean)=Top.index
				#-----------debug-----------#
				#print(round(cor.mean, 4))
				#---------------------------#
			}
			
			if(qtn.model == "SR"){
				if((is.null(inc.QTN.store.k) & is.null(inc.QTN.store.qtn))){
					model.index <- 1
				}else{
					if(mean(cor.store.k) > mean(cor.store.qtn)){
						model.index <- 2
					}else{
						model.index <- 3
					}
				}
			}
			
			model <- ifelse((is.null(inc.QTN.store.k) & is.null(inc.QTN.store.qtn))|model.index == 1, "K", ifelse(model.index == 2, "QTN+K", "QTN"))
			
			if(model == "K"){
				cross.QTN = NULL
				cross.model=model
			}
			if(model == "QTN+K"){
				cross.QTN=inc.QTN.store.k
				cross.model=model
			}
			if(model == "QTN"){
				cross.QTN=inc.QTN.store.qtn
				cross.model=model
			}
			
			cat("\r", "The Posterior QTNs:", paste(rep(" ", 20), collapse=""), "\n")
			if(is.null(cross.QTN)){
				cat(" NULL", "\n")
			}else{
				cat(" ", cross.QTN, "\n")
			}
			cat(" The optimized MODEL:", "\n")
			cat(" ", cross.model, "\n")
		}

	}
	#-----------debug-----------#
	# cpu=1; setMKLthreads(30)
	#---------------------------#
	
	sample.num <- sam
	
	if((!is.null(Top.perc) && (length(Top.perc) > 1 | length(amplify) > 1)) & cross.model!="QTN"){
		cat(" Optimizing KINSHIP...\n")
		if(K.method %in% c("center", "scale")){
			SUM <- n.marker
		}
		if(K.method == "vanraden"){
			geno.matrix <- as.matrix(geno)
			Pi <- 0.5 * rowMeans(geno.matrix)
			SUM <- sum(Pi * (1-Pi))
			rm(geno.matrix); rm(Pi); gc()
		}
		
		Top.perc <- sort(Top.perc)
		if(!0 %in% Top.perc)	Top.perc <- c(0, Top.perc)
		
		mult.run <- function(j, math.cpu=NULL){
			if(cpu>1 & (wind | r.open)) {
				try(setMKLthreads(math.cpu), silent=TRUE)
			}
			cv <- j %% (length(Top.perc) * length(amplify))
			if(cv == 0){
				cv <- j %/% (length(Top.perc) * length(amplify))
				top <- length(Top.perc)
				amp.index <- length(amplify)
			}else					
			{
				cv <- j %/% (length(Top.perc) * length(amplify)) + 1
				amp.index <- (j %% (length(Top.perc) * length(amplify))) %% length(amplify)
				if(amp.index == 0){
					amp.index <- length(amplify)
					top <- (j %% (length(Top.perc) * length(amplify))) %/% length(amplify)
				}else
				{
					top <- (j %% (length(Top.perc) * length(amplify))) %/% length(amplify) + 1
				}
			}
			
			mult.inf.index <- inf.index[[cv]]
			mult.ref.index <- !NA.index
			mult.ref.index[mult.inf.index] <- FALSE
			
			P.ref <- phe
			P.ref[!mult.ref.index] <- NA
			P.inf <- phe[mult.inf.index]

			if(Top.perc[top] == 0){
				if(amp.index == 1){
					Kt <- K
					max.wt <- 1
				}else
				{
					Kt <- NULL
				}
			}else
			{
				mytop.order <- order(P.value[[cv]], decreasing=FALSE)
				mytop.perc <- mytop.order[1:ceiling(n.marker * Top.perc[top])]
                mytop.sub <- mytop.order[ceiling(n.marker * Top.perc[top]) + 1]
                top.wt <- log(P.value[[cv]][mytop.sub], base=amplify[amp.index])-log(P.value[[cv]][mytop.perc], base=amplify[amp.index])
				# top.wt <- -log(P.value[[cv]][mytop.perc], base=amplify[amp.index])
				# top.wt <- -log(P.value[[cv]][mytop.perc])
				max.wt <- max(top.wt)
				
				#set the weight of cross.QTN to 0
				#top.wt[cross.QTN] <- 0
				
				#-----------debug-----------#
                #print(range(top.wt))
                #print(sum(is.na(top.wt)))
				#---------------------------#

				# Kt <- KAML.Kin(geno[mytop.perc, ], weight=top.wt, SUM=SUM, type=K.method)
				Kt <- KAML.KinNew(deepcopy(geno, rows=mytop.perc), SUM=SUM, scale=(K.method=="scale"), priority=priority, weight=top.wt, threads=1)
				Kt <- (K + Kt)/(mean(diag(K)) + mean(diag(Kt)))
				# Kt <- (1 - amplify[amp.index]) * K + (amplify[amp.index]) * Kt
				
				
				#-----------debug-----------#
				 # write.csv(P.value[[cv]], "debug.csv")
				  # print(sum(P.value[[cv]] == 0))
				  # print(P.value[[cv]][mytop.sub])
				  # print(range(P.value[[cv]][mytop.perc]))
				  # print(Kt[1:3,1:3])
				  # print(range(top.wt))
				# if(j <= 20)	write.csv(myK, paste(j,"_",amplify[amp.index], "_err.K.csv", sep=""), row.names=FALSE)
				#---------------------------#

				rm(list=c("mytop.order", "mytop.perc", "top.wt")); gc()
			}
			if(!is.null(Kt)){
				if(!is.null(cross.QTN)){
					QTN.cv <- t(geno[cross.QTN, , drop=FALSE])
				}else
				{
					QTN.cv <- NULL
				}
				KAML.cv <- cbind(CV, QTN.cv)
				gblup <- KAML.Mix(phe=P.ref, vc.method=vc.method, CV=KAML.cv, K=Kt)
				KAML.ll <- gblup$LL
				gebv <- (gblup$ebv + KAML.cv %*% as.matrix(gblup$beta))[mult.inf.index]
				if(binary){
					acc<- KAML.stas.cal(P.inf, gebv, type="auc")
				}else{
					acc<- KAML.stas.cal(P.inf, gebv, type="cor")
				}
				rm(list=c("P.ref", "P.inf", "Kt", "gblup", "gebv")); gc()
				#cat(paste(" Cross-validation NO.", cv, "; QTN=", ceiling(n.marker * Top.perc[top]), "(", Top.perc[top] * 100,
                #"%); Logx=", round(amplify[amp.index], 4),"; ",mytext , "=", round(acc, 4), paste(rep(" ", 5), collapse=""), "\r", sep=""))
			}else{
				acc <- NA; max.wt <- NA; KAML.ll <- NA
			}
			print.f(j)
			return(list(acc=acc,max.wt=max.wt,KAML.ll=KAML.ll))
		}
		iterationN <- ((sample.num * crv.num) * length(Top.perc) * length(amplify))
		cat(" Grid search started...\n")
		cat(paste(" Total iteration number:", iterationN))
		if(cpu == 1){
			print.f <- function(i){KAML.Bar(i=i, n=iterationN, type="type1", symbol.len=0, symbol.head=" Cross-validation Finished_", symbol.tail="")}
			mult.res <- lapply(1 : iterationN, mult.run)
		}else
		{
			if(wind){
				print.f <- function(i){KAML.Bar(i=i, n=iterationN, type="type1", symbol.len=0, symbol.head=" Cross-validation Finished_", symbol.tail="")}
				cat(" Multi-process started...\n")
				cat(" (Note: There needs to wait some time! See iteration details in 'Loop.log')\n")
				if(is.null(Top.index)){
					cl <- makeCluster(cpu, outfile = "Loop.log")
					registerDoParallel(cl)
					clusterExport(cl, varlist=c("KAML.Mix", "KAML.EMMA.REML", "KAML.stas.cal", "KAML.Kin", "KAML.QTN.rm", "KAML.GLM"))
					#Exp.packages <- clusterEvalQ(cl, c(library(bigmemory), library(rfunctions)))
				}
				mult.res <- foreach(x=1:iterationN,
                .packages=c("bigmemory", "rfunctions")) %dopar% mult.run(x, math.cpu=mkl)
				if(!bisection)	stopCluster(cl)
				cat(" Multi-process done!\n")
			}else
			{
				tmpf.name <- tempfile()
				tmpf <- fifo(tmpf.name, open="w+b", blocking=TRUE)
				
				# print.f <- function(i){writeBin(1, tmpf)}
				# KAML.Bar(n=iterationN, type="type2", tmp.file=tmpf, fixed.points=FALSE, symbol.len=0, symbol.head=" Cross-validation Finished_", symbol.tail="")
				
				writeBin(0, tmpf)
				print.f <- function(i){KAML.Bar(n=iterationN, type="type3", tmp.file=tmpf, fixed.points=FALSE, symbol.len=0, symbol.head=" Cross-validation Finished_", symbol.tail="")}
				
				if(r.open)	setMKLthreads(1)
				mult.res <- parallel::mclapply(1:iterationN, mult.run, mc.cores = cpu)
				if(r.open)	setMKLthreads(mkl)
				close(tmpf); unlink(tmpf.name)
			}
		}
		
		mult.store <- do.call(rbind, mult.res)
		K.cor.store <- matrix(unlist(mult.store[, 1]), length(Top.perc) * length(amplify))
		K.cor.store <- apply(K.cor.store, 2, function(x){x[which(is.na(x))]=x[1]; return(x)})
		Max.wt.store <- matrix(unlist(mult.store[, 2]), length(Top.perc) * length(amplify))
		LL.store <- matrix(unlist(mult.store[, 3]), length(Top.perc) * length(amplify))
		cor_qtn_k <- K.cor.store[1, ]
        rm(list=c("mult.store", "mult.res")); gc()

        top.index <- rep(Top.perc, rep(length(amplify), length(Top.perc)))
        top.amplify <- rep(amplify, length(Top.perc))
		colnames(K.cor.store) <- 1:(sample.num * crv.num)
		rownames(K.cor.store) <- paste(top.index, top.amplify,sep="_")
        rownames(Max.wt.store) <- paste(top.index, top.amplify,sep="_")
		rownames(LL.store) <- paste(top.index, top.amplify,sep="_")
		
        #-----------debug-----------#
        #write.csv(K.cor.store,"K.store.debug.csv")
        #print(K.cor.store)
        #print(Max.wt.store)
		#print(LL.store)
        #---------------------------#
		
		K.cor.store <- K.cor.store[-c(1:length(amplify)), ]
		Top.perc <- Top.perc[-1]
        top.index <- rep(Top.perc, rep(length(amplify), length(Top.perc)))
        top.amplify <- rep(amplify, length(Top.perc))
		colnames(K.cor.store) <- 1:(sample.num * crv.num)
		rownames(K.cor.store) <- paste(top.index, top.amplify,sep="_")
		
		K.cor.mean <- apply(K.cor.store, 1, mean)
		#K.cor.mean <- apply(K.cor.store, 1, median)
		top.index.max <- top.index[which.max(K.cor.mean)]
		amp.max <- top.amplify[which.max(K.cor.mean)]
		cat("\r", paste("Prior top percentage: ", paste(rep(" ", 20), collapse=""), "\n", sep=""))
		cat(" ", paste(top.index.max * 100, "%", sep=""), "\n")
		cat(paste(" Prior Logx: ", "\n", sep=""))
		cat(" ", round(amp.max, 4), "\n")
		
		bisection <- ifelse(bisection.loop != 0, TRUE, FALSE)
		if(bisection){
			cat(" Bisection algorithm started...\n")
			TOP <- Top.perc
			AMP <- amplify
			K.COR <- K.cor.store
			for(loop in 1 : bisection.loop){
				K.cor.mean <- apply(K.cor.store, 1, mean)
				top.cor <- K.cor.store[which.max(K.cor.mean), ]
				K.cor.max.mean <- max(K.cor.mean)
				if(length(Top.perc) != 1){
					top.index.max <- top.index[which.max(K.cor.mean)]
					max.pos <- which(Top.perc == top.index.max)
					#K.top.max <- tapply(K.cor.mean, top.index, max)
					#K.top.order <- order(K.top.max, decreasing=TRUE)
					top <- Top.perc[max.pos]
					TOP.index <- which(TOP == top)
					if(TOP.index == 1){
						flank1=0; flank2=TOP[TOP.index+1]
					}else if(TOP.index == length(TOP)){
						flank1=TOP[TOP.index-1]; flank2=2 * TOP[TOP.index] - TOP[TOP.index-1]
					}else{
						flank1=TOP[TOP.index-1]; flank2=TOP[TOP.index+1]
					}
					top1 <- mean(c(flank1, top))
					top2 <- mean(c(top, flank2))
					TOP <- c(top1, top2, TOP)
					TOP <- sort(TOP)
					Top.perc <- c(top1, top2)
				}else{
					top1 <- Top.perc
					top2 <- Top.perc
				}
				
				if(length(amplify) != 1){
					amp.max <- top.amplify[which.max(K.cor.mean)]
					max.pos <- which(amplify == amp.max)
					amp <- amplify[max.pos]
					AMP.index <- which(AMP == amp)
					if(AMP.index == 1){
						flank1=1; flank2=AMP[AMP.index+1]
					}else if(AMP.index == length(AMP)){
						flank1=AMP[AMP.index-1]; flank2=2 * AMP[AMP.index] - AMP[AMP.index-1]
					}else{
						flank1=AMP[AMP.index-1]; flank2=AMP[AMP.index+1]
					}
					amp1 <- mean(c(flank1, amp))
					amp2 <- mean(c(amp, flank2))
					AMP <- c(amp1, amp2, AMP)
					AMP <- sort(AMP)
					amplify <- c(amp1, amp2)
				}else{
					amp1 <- amplify
					amp2 <- amplify
				}
				iterationN <- ((sample.num * crv.num) * length(Top.perc) * length(amplify))
				cat("\r", paste("Loop ", loop, " of ", bisection.loop, " selected interval: [", round(top1 * 100, 4), "%, ", round(top2 * 100, 4), "%]; ", "[", round(amp1, 4), ", ", round(amp2, 4), "]\n", sep=""))
				cat(paste(" Total iteration number:", iterationN))
				if(cpu == 1){
					print.f <- function(i){KAML.Bar(i=i, n=iterationN, type="type1", symbol.len=0, symbol.head=" Cross-validation Finished_", symbol.tail="")}
					mult.res <- lapply(1:iterationN, mult.run)
				}else
				{
					if(wind){
						print.f <- function(i){KAML.Bar(i=i, n=iterationN, type="type1", symbol.len=0, symbol.head=" Cross-validation Finished_", symbol.tail="")}
						mult.res <- foreach(x=1:iterationN,
						.packages=c("bigmemory", "rfunctions")) %dopar% mult.run(x, math.cpu=mkl)
						if(loop == bisection.loop)	stopCluster(cl)
					}else
					{
						tmpf.name <- tempfile()
						tmpf <- fifo(tmpf.name, open="w+b", blocking=TRUE)
						
						# print.f <- function(i){writeBin(1, tmpf)}
						# KAML.Bar(n=iterationN, type="type2", tmp.file=tmpf, fixed.points=FALSE, symbol.len=0, symbol.head=" Cross-validation Finished_", symbol.tail="")
						
						writeBin(0, tmpf)
						print.f <- function(i){KAML.Bar(n=iterationN, type="type3", tmp.file=tmpf, fixed.points=FALSE, symbol.len=0, symbol.head=" Cross-validation Finished_", symbol.tail="")}
				
						if(r.open)	setMKLthreads(1)
						mult.res <- parallel::mclapply(1:iterationN, mult.run, mc.cores = cpu)
						if(r.open)	setMKLthreads(mkl)
						close(tmpf); unlink(tmpf.name)
					}
				}
				mult.store <- do.call(rbind, mult.res)
				K.cor.store <- matrix(unlist(mult.store[, 1]), length(Top.perc) * length(amplify))
				K.COR <- rbind(K.COR, K.cor.store)
				K.cor.store <- rbind(top.cor, K.cor.store)
				if(length(amplify) == 1){
					top.index <- c(top, Top.perc)
				}else if(length(Top.perc) == 1){
					top.index <- c(top, rep(Top.perc, 2))}else{top.index <- c(top, rep(Top.perc, c(2, 2)))
				}
				top.amplify <- c(amp, rep(amplify, length(Top.perc)))
				colnames(K.cor.store) <- 1:(sample.num * crv.num)
				rownames(K.cor.store) <- paste(top.index, top.amplify,sep="_")
				if(length(Top.perc) != 1)	Top.perc <- c(top, Top.perc)
				if(length(amplify) != 1)	amplify <- c(amp, amplify)
				
				#-----------debug-----------#
				# print(K.cor.store)
				#---------------------------#
				
				rm(list=c("mult.store", "mult.res")); gc()
				
				#test converge
				K.cor.mean <- apply(K.cor.store, 1, mean)
				top.index.max <- top.index[which.max(K.cor.mean)]
				amp.max <- top.amplify[which.max(K.cor.mean)]
				diff <- (max(K.cor.mean) - K.cor.max.mean)
				if(diff > 0 & diff < 1e-5){
					cat("\r Loop Stopped at convergence!                          \n")
					break()
				}
			}
		}
		
        if(
			(max(K.cor.mean) > mean(cor_qtn_k))
			#& 
			#sum((K.cor.store[which.max(K.cor.mean), ] - cor_qtn_k) > 0) >= floor((sample.num * crv.num) * count.threshold)
		)
		{
			K.opt <- TRUE
			cat("\r", paste("Posterior top percentage: ", paste(rep(" ", 20), collapse=""), "\n", sep=""))
			cat(" ", paste(top.index.max * 100, "%", sep=""), "\n")
			cat(paste(" Posterior Logx: ", "\n", sep=""))
			cat(" ", round(amp.max, 4), "\n")
		}else{
			K.opt <- FALSE
			top.index.max <- NULL
			amp.max <- NULL
			cat("\r", paste("Posterior top percentage:", paste(rep(" ", 20), collapse=""), "\n", sep=""))
			cat(" NULL", "\n")
			cat(paste(" Posterior Logx:", "\n", sep=""))
			cat(" NULL", "\n")
		}

		if(r.open & linux & cpu > 1)	try(setMKLthreads(mkl), silent=TRUE)
		
		if(K.opt){
			mytop.order <- order(P.value.ref, decreasing = FALSE)
			mytop.num <- ceiling(n.marker * top.index.max)
			mytop.perc <- mytop.order[1 : mytop.num]
            mytop.sub <- mytop.order[mytop.num + 1]
            
            top.wt <- log(P.value.ref[mytop.sub], base=amp.max) - log(P.value.ref[mytop.perc], base=amp.max)
			# top.wt <- -log(P.value.ref[mytop.perc], base=amp.max)
			# top.wt <- -log(P.value.ref[mytop.perc])
			
			#top.wt <- top.wt*1.5
			
			#-----------debug-----------#
            # print(range(top.wt))
			#---------------------------#
			
			cat(" Calculating optimized KINSHIP...\n")
			# optK <- KAML.Kin(geno[mytop.perc,], type=K.method, weight = top.wt, SUM=SUM)
            optK <- KAML.KinNew(deepcopy(geno, rows=mytop.perc), SUM=SUM, scale=(K.method=="scale"), priority=priority, weight=top.wt, threads=cpu, verbose=TRUE)
				
            K <- (K + optK)/(mean(diag(K)) + mean(diag(optK)))
            # Kt <- (1 - amp.max) * K + (amp.max) * optK
				
            
            #-----------debug-----------#
            #print(K[1:5,1:5])
            #---------------------------#
            
			rm(list=c("mytop.order", "mytop.perc", "top.wt", "optK"))
		}
		rm(list=c("P.value", "P.value.ref")); gc()
	}else
	{
		if((length(Top.perc) == 1 & length(amplify) == 1) & cross.model!="QTN"){
			cat(" The provided Top percentage: \n")
			cat(" ", paste(Top.perc * 100, "%\n", sep=""))
			cat(" The provided Logx: \n")
			cat(" ", amplify, "\n")
			top.index.max <- Top.perc
			amp.max <- amplify
			
			mytop.order <- order(P.value.ref, decreasing = FALSE)
			mytop.num <- ceiling(n.marker * Top.perc)
			mytop.perc <- mytop.order[1 : mytop.num]
            mytop.sub <- mytop.order[mytop.num + 1]
            
            top.wt <- log(P.value.ref[mytop.sub], base=amplify) - log(P.value.ref[mytop.perc], base=amplify)
			
			#top.wt <- top.wt*1.5
			
			#-----------debug-----------#
            # print(range(top.wt))
			#---------------------------#
			
			cat(" Calculating optimized KINSHIP...\n")
			# optK <- KAML.Kin(geno[mytop.perc,], type="center", weight = top.wt, SUM=nrow(geno))
            optK <- KAML.KinNew(deepcopy(geno, rows=mytop.perc), SUM=nrow(geno), scale=(K.method=="scale"), priority=priority, weight=top.wt, threads=cpu, verbose=TRUE)
			
            K <- (K + optK)/(mean(diag(K)) + mean(diag(optK)))
			rm(list=c("mytop.order", "mytop.perc", "top.wt", "optK"))
			rm(list=c("P.value.ref")); gc()
		}else{
			cat(" No Kinship optimization\n")
			top.index.max <- NULL
			amp.max <- NULL
		}
	}
	cat(" Cross-validation DONE!\n")
	return(list(cross.QTN = cross.QTN, cross.model=cross.model, cross.tp=top.index.max, cross.amp=amp.max, cross.k=K))
}

KAML.GLM <- 
function(
	y, X=NULL, qtn.matrix
)
{
#--------------------------------------------------------------------------------------------------------#
# Object: To solve the GLM model and get the effects of QTNs											 #
#																										 #
# Input:	 																							 #
# y: a vector of phenotype	 																			 #
# X: The fixed effect(X must contain a column of 1's)													 #
# qtn.matrix: a  n1 * m1 matrix, n1 is the population size, m1 is the number of selected QTNs			 #
#--------------------------------------------------------------------------------------------------------#
	r.open <- !inherits(try(Revo.version,silent=TRUE),"try-error")
	qtn.matrix <- data.matrix(qtn.matrix)
	if(!is.numeric(y))	y <- as.numeric(as.character(y))
	Y<-matrix(y)
	X<-cbind(X, qtn.matrix)
	XX <- KAML.Crossprod(X)
	XY<-crossprod(X, Y)
	YY <- KAML.Crossprod(Y)
	XXi <- try(solve(XX), silent = TRUE)
	if(inherits(XXi, "try-error")){
		XXi <- KAML.ginv(XX)
	}
	beta<-crossprod(XXi, XY)
	QTN.eff<-as.numeric(beta)
	return(list(QTN.eff=QTN.eff))
}

KAML.GWAS <- 
function(
	phe, geno, K=NULL, CV=NULL, NPC=NULL, REML=NULL, cpu=1, priority="speed", vc.method="emma", method="MLM", bar.head="|", bar.tail=">", bar.len=50
)
{
#--------------------------------------------------------------------------------------------------------#
# Object: To perform GWAS using GLM or MLM model, get the P value of SNPs								 #
#																										 #
# Input:																								 #
# phe: Phenotype, a value vector(NA is allowed, only non-NA individuals will be used for analysis)		 #
# geno: Genotype in numeric format, pure 0, 1, 2 matrix; m * n, m is marker size, n is population size	 #
# (Note: both matrix or big.matrix are acceptable)														 #
# K: Kinship for all individuals(row and column orders must be the same with phenotype)					 #
# CV: covariance, design matrix(n * x) for the fixed effects(CV must contain a column of 1's)			 # 
# NPC: the number of PC that will be added as covariance to control population structure				 #
# REML: a list contains ve and vg																		 #
# vc.method: method for variance components estimation("emma" or "gemma")								 #
# cpu: the number of CPU for calculation																 #
# method: "GLM" or "MLM"	 																			 #
# bar.head: the head of Progress bar	 																 #
# bar.len: the length of the bar	 																	 #
#--------------------------------------------------------------------------------------------------------#
	R.ver <- Sys.info()[['sysname']]
	wind <- R.ver == 'Windows'
	linux <- R.ver == 'Linux'
	mac <- (!linux) & (!wind)
    r.open <- !inherits(try(Revo.version,silent=TRUE), "try-error")
	math.cpu <- Math_cpu_check()
	if(r.open &  mac & cpu > 1)	Sys.setenv("VECLIB_MAXIMUM_THREADS" = "1")
	if(wind)	cpu <- 1

	bar <- ifelse(bar.len == 0, FALSE, TRUE)
	
	Ind.index <- !is.na(phe)
	ys <- as.numeric(phe[Ind.index])
	if(sum(Ind.index) == length(phe)){
		# geno <- as.matrix(geno)
	}else{
		geno <- deepcopy(geno, cols=Ind.index)
	}
	n <- ncol(geno)
	m <- nrow(geno)
	if(method == "MLM"){
		if(is.null(K)){
			# K <- KAML.Kin(geno, priority=priority)
			K <- KAML.KinNew(geno, priority=priority, threads=cpu)
			
		}else{
			K <- K[Ind.index, Ind.index]
		}
	}
	if(method == "MLM"){
		mkl.cpu <- ifelse((2^(n %/% 1000)) < math.cpu, 2^(n %/% 1000), math.cpu)
		try(setMKLthreads(mkl.cpu), silent=TRUE)
        eig <- eigen(K, symmetric=TRUE)
		try(setMKLthreads(math.cpu), silent=TRUE)
	}else{
		eig <- NULL
	}
	if(!is.null(NPC)){
		if(!is.null(eig)){
			pca <- eig$vectors[, c(1:NPC)]
		}else{
			pca <- prcomp((t(geno)))$x[, 1:NPC]
		}
	}else{
		pca <- NULL
	}
	if(is.null(CV)){
		X0 <- cbind(matrix(1, n), pca)
	}else{
		X0 <- cbind(matrix(1, n), pca, CV[Ind.index, ])
	}
	if(is.null(REML) & method == "MLM"){
		if(vc.method == "he") REML <- KAML.HE(ys, X=X0, K=K)
		if(vc.method == "emma") REML <- KAML.EMMA.REML(ys, X=X0, K=K)
		if(vc.method == "brent") REML <- KAML.EIGEN.REML(ys, X=X0, eigenK=eig)
	}
	# q0 <- ncol(X0)
	# iXX <- matrix(NA,q0+1,q0+1)
	# Xt <- matrix(NA,n, q0+1)
	
	#parallel function for MLM model
	eff.mlm <- function(i){
		SNP <- as.matrix(geno[i, ])
		xst <- crossprod(U, SNP)
		Xt[1:n,q0+1] <- xst
		X0Xst <- crossprod(X0t,xst)
		XstX0 <- t(X0Xst)
		xstxst <- crossprod(xst, xst)
		xsY <- crossprod(xst,yt)
		XY <- c(X0Y,xsY)
		B22 <- xstxst - XstX0%*%iX0X0%*%X0Xst
		invB22 <- 1/B22
		B21 <- tcrossprod(XstX0, iX0X0)
		NeginvB22B21 <- crossprod(-invB22,B21)
		B11 <- iX0X0 + as.numeric(invB22)*crossprod(B21,B21)
		iXX[1:q0,1:q0]=B11
		iXX[(q0+1),(q0+1)]=1/B22
		iXX[(q0+1),1:q0]=NeginvB22B21
		iXX[1:q0,(q0+1)]=NeginvB22B21
		beta <- crossprod(iXX,XY)
		stats <- beta[(q0+1)]/sqrt((iXX[(q0+1), (q0+1)]) * vgs)
		p <- 2 * pt(abs(stats), n-(q0+1), lower.tail=FALSE)
		effect<- beta[(q0+1)]
		if(bar)	print.f(i)
		return(list(effect = effect, p = p))
	}
	
	#parallel function for GLM model
	eff.glm <- function(i){
		SNP <- geno[i, ]
		sy <- crossprod(SNP,y)
		ss <- crossprodcpp(SNP)
		xs <- crossprod(X0,SNP)
		B21 <- crossprod(xs, X0X0i)
		t2 <- B21 %*% xs
		B22 <- ss - t2
		invB22 <- 1/B22
		NeginvB22B21 <- crossprod(-invB22,B21)
		B21B21 <- crossprodcpp(B21)
		iXX11 <- X0X0i + as.numeric(invB22) * B21B21
		iXX[1:q0,1:q0] <- iXX11
		iXX[(q0+1),(q0+1)] <- invB22
		iXX[(q0+1),1:q0] <- NeginvB22B21
		iXX[1:q0,(q0+1)] <- NeginvB22B21
		df <- n-q0-1
		rhs <- c(X0Y,sy)
		effect <- crossprod(iXX,rhs)
		ve <- (YY-crossprod(effect,rhs))/df
		effect <- effect[q0+1]
		t.value <- effect/sqrt(iXX[q0+1, q0+1] * ve)
		p <- 2 * pt(abs(t.value), df, lower.tail=FALSE)
		if(bar)	print.f(i)
		return(list(effect=effect, p=p))
    }
    
    if(method == "MLM"){
        ves <- REML$ve
        vgs <- REML$vg
        lambda <- ves/vgs
        U <- eig$vectors * matrix(sqrt(1/(eig$values + lambda)), n, length(eig$values), byrow=TRUE); rm(eig); gc()
		# y <- matrix(ys)
		# yt <- crossprod(U, y)
		# X0t <- crossprod(U, X0)
		# X0X0 <- KAML.Crossprod(X0t)
		# X0Y <- crossprod(X0t, yt)
		# iX0X0 <- ginv(X0X0)
		# Xt[1:n,1:q0] <- X0t
    }
	
  #   if(method == "GLM"){
		# y <- matrix(ys)
  #       X0X0 <- KAML.Crossprod(X0)
		# X0Y <- crossprod(X0, y)
		# YY <- KAML.Crossprod(y)
		# X0X0i <- ginv(X0X0)
  #   }
	
	# options(warn = -1)
 #    if(cpu == 1){
	# 	if(bar)	print.f <- function(i){KAML.Bar(i=i, n=m, type="type1", symbol.len=bar.len, symbol.head=bar.head, symbol.tail=bar.tail)}
	# 	mkl.cpu <- ifelse((2^(n %/% 1000)) < math.cpu, 2^(n %/% 1000), math.cpu)
	# 	if(r.open)	try(setMKLthreads(mkl.cpu), silent=TRUE)
	# 	if(method == "MLM") results <- lapply(1:m, eff.mlm)
	# 	if(method == "GLM")	results <- lapply(1:m, eff.glm)
	# 	if(r.open)	try(setMKLthreads(math.cpu), silent=TRUE)
	# }else
	# {
	# 	if(method == "MLM"){
	# 		if(wind){
	# 			if(bar)	print.f <- function(i){KAML.Bar(i=i, n=m, type="type1", symbol.len=bar.len, symbol.head=bar.head, symbol.tail=bar.tail)}
	# 			#foreach function
	# 			#print(" *  *  *  *  *  * test foreach parallel *  *  *  *  *  * ")
	# 			#cl <- makeCluster(cpu)
	# 			#registerDoParallel(cl)
	# 			#Exp.packages <- clusterEvalQ(cl, c(library(bigmemory)))
	# 			#results <- foreach(x=1:m) %dopar% eff.mlm(x)
	# 			#stopCluster(cl)
				
	# 			#print(" *  *  *  *  *  * test parLapply parallel *  *  *  *  *  * ")
	# 			cl <- makeCluster(getOption("cl.cores", cpu))
	# 			clusterExport(cl, varlist=c("geno", "yt", "X0", "U", "vgs", "ves", "mkl"), envir=environment())
	# 			Exp.packages <- clusterEvalQ(cl, c(library(bigmemory),library(rfunctions)))
	# 			results <- parallel::parLapply(cl, 1:m, eff.mlm)
	# 			stopCluster(cl)
	# 		}else{
	# 			if(bar){
	# 				tmpf.name <- tempfile()
	# 				tmpf <- fifo(tmpf.name, open="w+b", blocking=TRUE)
	# 				writeBin(0, tmpf)
	# 				print.f <- function(i){KAML.Bar(n=m, type="type3", tmp.file=tmpf, symbol.len=bar.len, symbol.head=bar.head, symbol.tail=bar.tail)}		
	# 			}
	# 			if(r.open){
	# 				if(R.ver == 'Linux' && r.open)	try(setMKLthreads(1), silent=TRUE)
	# 				results <- parallel::mclapply(1:m, eff.mlm, mc.cores=cpu)
	# 				if(R.ver == 'Linux' && r.open)	try(setMKLthreads(math.cpu), silent=TRUE)
	# 			}else{
	# 				results <- parallel::mclapply(1:m, eff.mlm, mc.cores=cpu)
	# 			}
	# 			#Sys.sleep(0.5)
	# 			if(bar){close(tmpf); unlink(tmpf.name); cat('\n');}
	# 		}
	# 	}
	# 	if(method  == "GLM"){
	# 		if(wind){
	# 			if(bar)	print.f <- function(i){KAML.Bar(i=i, n=m, type="type1", symbol.len=bar.len, symbol.head=bar.head, symbol.tail=bar.tail)}
	# 			#parLapply function
	# 			#print(" *  *  *  *  *  * test parLapply parallel *  *  *  *  *  * ")
	# 			cl <- makeCluster(getOption("cl.cores", cpu))
	# 			clusterExport(cl, varlist=c("geno", "ys", "X0", "math.cpu"), envir=environment())
	# 			Exp.packages <- clusterEvalQ(cl, c(library(bigmemory),library(rfunctions)))
	# 			results <- parallel::parLapply(cl, 1:m, eff.glm)
	# 			stopCluster(cl)
	# 		}else{
	# 			if(bar){
	# 				tmpf.name <- tempfile()
	# 				tmpf <- fifo(tmpf.name, open="w+b", blocking=T)
	# 				writeBin(0, tmpf)
	# 				print.f <- function(i){KAML.Bar(n=m, type="type3", tmp.file=tmpf, symbol.len=bar.len, symbol.head=bar.head, symbol.tail=bar.tail)}
	# 			}
	# 			if(r.open){
	# 				if(R.ver == 'Linux' && r.open)	try(setMKLthreads(1), silent=TRUE)
	# 				results <- parallel::mclapply(1:m, eff.glm, mc.cores=cpu)
	# 				if(R.ver == 'Linux' && r.open)	try(setMKLthreads(math.cpu), silent=TRUE)
	# 			}else{
	# 				results <- parallel::mclapply(1:m, eff.glm, mc.cores=cpu)
	# 			}
	# 			#Sys.sleep(0.5)
	# 			if(bar){close(tmpf); unlink(tmpf.name); cat('\n');}
	# 		}
	# 	}
	# }
	# options(warn = 0)

	if(method == "GLM"){
		if(cpu > 1 & r.open)	try(setMKLthreads(1), silent=TRUE)
		results <- glm_c(y=ys, X=X0, geno@address, barhead=bar.head, verbose=bar, threads=cpu)
		if(cpu > 1 & r.open)	try(setMKLthreads(math.cpu), silent=TRUE)
	}else if(method == "MLM"){
		if(cpu > 1 & r.open)	try(setMKLthreads(1), silent=TRUE)
		results <- mlm_c(y=ys, X=X0, U=U, vgs=vgs, geno@address, barhead=bar.head, verbose=bar, threads=cpu)
		if(cpu > 1 & r.open)	try(setMKLthreads(math.cpu), silent=TRUE)
	}else{
		stop("Unknown gwas model.")
	}

	return(results)
}

KAML.copy <- 
function(
	x, cols=NULL, rows=NULL, memo="", backed=TRUE, delete=FALSE
)
{
#--------------------------------------------------------------------------------------------------------#
# Object: copy big.matrix and write the file to disk													 #
#																										 #
# Input:																								 #
# x: the big matrix which needs to be copied															 #
# cols: selected columns																				 #
# rows: selected rows																					 #
# memo: the specific names for output files																 #
# backed: whether to write the files to disk															 #
# delete: whether to delete the writed files in disk													 #
#--------------------------------------------------------------------------------------------------------#
	if(!is.big.matrix(x)) stop("x must be big.matrix!")
	if(backed){
	
		if(memo == ""){
			memo <- "KAML.temp"
		}else{
			memo <- paste(memo, ".temp", sep="")
		}

		copy.r <- deepcopy(x, cols=cols, rows=rows, backingfile=paste(memo, ".bin", sep=""), descriptorfile=paste(memo, ".desc", sep=""))

		if(delete){
			unlink(c(paste(memo, ".bin", sep=""), paste(memo, ".desc", sep="")), recursive = TRUE)
		}
	}else{
		copy.r <- deepcopy(x, cols=cols, rows=rows)
	}
	return(copy.r)
}

KAML.AIEM <- 
function(
	y, X=NULL, R=NULL, K=NULL, eigenK=NULL, start=NULL, method=c("AI", "EM", "EMAI", "AIEM"), nEMiter=1, nAIiter=20, tol=1e-6, verbose=TRUE, cpu=NULL
)
{

#To: variance component estimation by "AIREML" and "EMREML"
#Author: Lilin Yin
#Time: 2018/07/18
#y: numeric vector
#X: fixed effect
#R: random effect
#K: list or matrix, single or multiple Kinships of genotyped individuals
#eigenK: list, the eigen of Kinship
#method: reml method
#nEMiter: the max number of EM step
#nAIiter: the max number of AI step
#tol: the stop standard
#verbose: wether to print iteration step
#cpu: the cpu number with MKL 

#NOTE: NAs are allowed in y, the order of individuals must be same between y and K
	
	if(inherits(try(Revo.version,silent=TRUE),"try-error"))	cpu=1
	if(!inherits(try(Revo.version,silent=TRUE),"try-error") && !is.null(cpu)){
		try(setMKLthreads(cpu),silent = TRUE)
	}
	
	nan_index <- !is.na(y)
	#sample size
	n <- sum(nan_index)
	Y <- as.matrix(y[nan_index])
	
	if(!is.null(K)){
	
		if(is.data.frame(K))	K <- as.matrix(K)
		if(is.matrix(K))	K <- list(K)
		lencheck <- lapply(K, function(x){if(nrow(x) != length(y) && ncol(x) != length(y))	stop("Rank of K not equals to length of y")})
		
		#scale diag of all Kinship to 1
		if(is.list(K))	K <- lapply(K, 
			function(x){
				if(!is.matrix(x)){
					x <- as.matrix(x)
					x <- x/mean(diag(x))
				}else{
					x <- x/mean(diag(x))
				}
			}
		)
		eigenK <- NULL
	}else{
		if(is.null(eigenK)){
			if(is.null(R))	stop("K or eigen of K must be provided!")
		}else{
			if(n != length(y))	stop("'NAs' is not allowed in y when eigenK provided!")
			if(!all(names(eigenK) == c("values", "vectors")))	stop("Not supporting multiple random when eigenK provided!")
			if(nrow(eigenK$vectors) != length(y))	stop("Rank of eigenK not equals to length of y")
			m <- 1
		}
	}
	
	#fixed effect
	if(!is.null(X)){
		X <- as.matrix(X)
		if(nrow(X) != length(y))	stop("Rows not match between y and X!")
		if(sum(is.na(X)) != 0)	stop("NAs are not allowed in X!")
		X.index <- apply(X, 2, function(x) length(unique(x)) > 1)
		X <- X[, X.index, drop=FALSE]
		X <- cbind(1, X)
		X <- X[nan_index, ]
	}else{
		X<- matrix(1, nrow(Y))
	}
	
	#random effect
	if(!is.null(R)){
		R.cov <- function(y){
			indexM <- data.frame(1:length(y), y)
			k <- diag(1, length(y))
			index <- t(
				do.call(
					cbind,tapply(indexM[, 1], indexM[, 2], function(x){
						if(length(x)>1){
							return(combn(x, 2))
						}else{
							return(NULL)
						}
					}
					)
				)
			)
			k[index] <- 1
			k[index[, c(2,1)]] <- 1
			return(k)
		}
		R <- as.matrix(R)
		if(nrow(R) != length(y))	stop("Rows not match between y and R!")
		if(sum(is.na(R)) != 0)	stop("NAs are not allowed in R!")
		if(!all(apply(R, 2, function(x) length(unique(x)) > 1)))	stop("Group number must be more than 1!")
		if(!all(apply(R, 2, function(x) length(unique(x)) < length(x))))	stop("Group number must be less than the row number of R!")
		R.k <- list()
		for(r in 1:ncol(R)){
			R.k[[r]] <- R.cov(R[, r, drop=TRUE])
		}
		K <- c(R.k, K)
	}else{
		R.k <- NULL
	}
	m <- length(K)
	
	# check for starting values
	if(!is.null(start)){
		if(length(start) != (m+1)){
			stop("length of start must equal to ncol(R) + length(K) + 1")
		}
	}
	
	switch(
		match.arg(method),
		"AI"={
			IterM <- rep("AI", nAIiter)
		},
		"EM"={
			IterM <- rep("EM", nEMiter)
		},
		"AIEM"={
			IterM <- rep(c("AI", "EM"), c(nAIiter, nEMiter))
		},
		"EMAI"={
			IterM <- rep(c("EM", "AI"), c(nEMiter, nAIiter))
		}
	)
	
	#return the prediction, so store original K
	if(!is.null(K)){
		K_ori <- K
		if(n != length(y)){
			for(i in 1:m){
				# subset matrix
				K_ori[[i]] <- K_ori[[i]][,nan_index]
				K[[i]] <- K[[i]][nan_index,nan_index]
			}
		}
		
		if(!is.null(cpu) && cpu == 1 && is.null(eigenK) && m == 1){
			if(verbose)	message("Vinv with Eigen decomposition!")
			for(i in 1:m){
				eigenK <- eigen(K[[i]], symmetric=TRUE)
			}
		}
		
	}else{
		if(verbose)	message("Vinv with Eigen decomposition!")
		K <- K_ori <- list(Vinv <- eigenK$vectors %*% diag(eigenK$values) %*% t(eigenK$vectors))
	}
	
	# initial values
	sigma2.p <- var(Y)
	tol <- as.vector(tol*sigma2.p)  # set convergence tolerance dependent on trait
	
	if(is.null(start)){
		sigma2.k <- rep((1/(m+1))*sigma2.p, (m+1))
	}else{
		sigma2.k <- as.vector(start)
	}
	zeroFLAG <- rep(FALSE, m+1) # indicator of elements that have converged to "0"
	
	sigma2.kplus1 <- NULL
	reps <- 0
	
	if(verbose){
		if(is.null(R)){
			cat("        ",
				paste(paste("Var_K", 1:length(K), "(SE)", sep="", collapse="       "), "Var_e(SE)",
				paste("h2_K", 1:length(K), "(SE)", sep="", collapse="       "), sep="       "),"\n", sep=""
			)
		}else{
			if(m == ncol(R)){
				cat("        ",
					paste(paste("Var_R",  1:ncol(R), "(SE)", sep="", collapse="       "), "Var_e(SE)",
					paste("h2_R", 1:ncol(R), "(SE)", sep="", collapse="       "), sep="       "),"\n", sep=""
				)
			}else{
				cat("        ",
					paste(paste(paste("Var_R", 1:ncol(R), "(SE)", sep="", collapse="       "), "       ", "Var_K", 1:(length(K)-ncol(R)), "(SE)", sep="", collapse="       "), "Var_e(SE)",
					paste(paste("h2_R", 1:ncol(R), "(SE)", sep="", collapse="       "), "       ", "h2_K", 1:(length(K)-ncol(R)), "(SE)", sep="", collapse="       "), sep="       "),"\n", sep=""
				)
			}
		}
	}
			
	repeat({
		reps <- reps+1

		# variance matrix
		if(is.null(eigenK)){
			V <- 0
			for(i in 1:m){
				V <- V + K[[i]]*sigma2.k[i]
			}
			diag(V) <- diag(V) + sigma2.k[m+1]
			
			# inverse
			#Vinv <- chol2inv(chol(V))
			Vinv <- solve(V)
			
		}else{
			
			#Vinv=U(1/(vg*D+ve))U
			eigenVinv <- tcrossprod(eigenK$vectors, diag(1/(sigma2.k[1] * eigenK$values + sigma2.k[2])))
			Vinv <- tcrossprod(eigenVinv, eigenK$vectors); rm(eigenVinv)

		}
		
		# projection matrix
		tXVinv <- crossprod(X, Vinv)
		XVinvXinv <- ginv(tXVinv %*% X)
		P <- Vinv - Vinv %*% X %*% XVinvXinv %*% tXVinv

		# matrices for later use
		PY <- crossprod(P,Y)
		PPY <- crossprod(P,PY)
		
		if(IterM[reps] == "AI"){
		
			# Average Information and Scores
			AI <- matrix(NA, nrow=(m+1), ncol=(m+1))
			score <- rep(NA,(m+1))        
			for(i in 1:m){
				PAPY <- crossprod(P,crossprod(K[[i]],PY))
				score[i] <- -0.5*(sum(P*K[[i]]) - crossprod(Y, PAPY)) 
				AI[i,i] <- 0.5*crossprod(PY, crossprod(K[[i]],PAPY)) # YPAPAPY
					if((i+1) <= m){
						for(j in (i+1):m){
							AI[i,j] <- 0.5*crossprod(PY, crossprod(K[[j]],PAPY)) # YPDPAPY, YPSPAPY
							AI[j,i] <- AI[i,j]
						}
					}
				AI[i,(m+1)] <- 0.5*crossprod(PY, crossprod(K[[i]],PPY)) # YPAPPY
				AI[(m+1),i] <- AI[i,(m+1)]
			}
			score[m+1] <- -0.5*(sum(diag(P)) - crossprod(Y, PPY))
			AI[(m+1),(m+1)] <- 0.5*crossprod(PY,PPY) # YPPPY        	
			
			# update
			AIinv <- try(solve(AI), silent=TRUE)
			if(inherits(AIinv, "try-error"))	AIinv <- MASS::ginv(AI)
			
			AIinvScore <- AIinv %*% score
			sigma2.kplus1 <- sigma2.k + AIinvScore
			
			sigma2.kplus1[zeroFLAG & sigma2.kplus1 < tol] <- 0 # set elements that were previously "0" and are still < 0 back to 0 (prevents step-halving due to this component)
			
			# step-halving if step too far
			k <- 1
			while(!all(sigma2.kplus1 >= 0)){
				k <- 0.5*k
				sigma2.kplus1 <- sigma2.k + k*AIinvScore
				sigma2.kplus1[zeroFLAG & sigma2.kplus1 < tol] <- 0
			}
			SE.cal <- function(vc, Vvc){
				VR <- function(v1, vv1, v2, vv2, c12){
					nv <- v2^2 * vv1 + v1^2 * vv2 - 2 * v1 * v2 * c12
					dv <- v2^4
					nv/dv
				}
				n <- length(vc)
				VarV <- diag(Vvc)
				P <- sum(vc)
				VarP <- sum(Vvc)
				GGE <- vc[-n]
				h2 <- GGE/P
				SEGE <- sqrt(VarV)
				SEh2 <- vector()
				for (j in 1:(n - 1)){
					SEh2[j] <- VR(GGE[j], VarV[j], P, VarP, sum(Vvc[j, ]))
				}
				SEh2 <- sqrt(SEh2)
				return(list(SEGE=SEGE, h2=h2, SEh2=SEh2))
			}
			SE.res <- SE.cal(sigma2.kplus1, AIinv)
			if(verbose){
				cat("[", IterM[reps], "] ",
					paste(paste(sprintf("%.6f", sigma2.kplus1), "(", sprintf("%.4f", SE.res$SEGE), ")", sep="", collapse=" "),
					paste(sprintf("%.4f", SE.res$h2), "(", sprintf("%.4f", SE.res$SEh2), ")", sep="", collapse=" "), sep=" "),"\n", sep=""
				)
			}
		}
		
		if(IterM[reps] == "EM"){
			for(i in 1:m){
				PAPY <- crossprod(P,crossprod(K[[i]],PY))
				sigma2.kplus1[i] <- (1/n)*((sigma2.k[i])^2*crossprod(Y,PAPY) + (n*sigma2.k[i] - (sigma2.k[i])^2*sum(P*K[[i]])))
			}
			sigma2.kplus1[m+1] <- (1/n)*((sigma2.k[m+1])^2*crossprod(Y,PPY) + (n*sigma2.k[m+1] - (sigma2.k[m+1])^2*sum(diag(P))))

			# print current estimates
			if(verbose) message("[", IterM[reps], "] ", paste(round(sigma2.kplus1, 4), collapse=" "))
		}

		# test for convergence
		#stat <- max(abs(sigma2.kplus1 - sigma2.k))
		stat <- sqrt(sum((sigma2.kplus1 - sigma2.k)^2))
		sigma2.k <- sigma2.kplus1
		if(stat < tol | reps == length(IterM)) break()
		zeroFLAG <- sigma2.k < tol # which elements have converged to "0"
		sigma2.k[zeroFLAG] <- 0 # set these to 0
	})

	beta <- XVinvXinv %*% tXVinv %*% Y
	u <- NULL
	for(i in 1:m){
		u <-  cbind(u, K_ori[[i]] %*% Vinv %*% (Y-X%*%beta) * sigma2.k[[i]])
	}
	
	return(list(vc = sigma2.k, h2=SE.res$h2, vcse=SE.res$SEGE, h2se=SE.res$SEh2, beta=beta, u=u))
}

KAML <- 
function(
	pfile="", gfile="", kfile=NULL, dcovfile=NULL, qcovfile=NULL, pheno=1, SNP.weight=NULL, GWAS.model=c("MLM","GLM", "RR"), GWAS.npc=NULL, prior.QTN=NULL, prior.model=c("QTN+K", "QTN", "K"), vc.method=c("brent", "he", "emma"), 
	K.method=c("center", "scale", "vanraden"), Top.perc=c(1e-4, 1e-3, 1e-2, 1e-1), Top.num=15, Logx=c(1.01, 1.11, exp(1), 10), qtn.model=c("MR", "SR", "BF"), BF.threshold=NULL, binary=FALSE,
	bin.size=1000000, max.nQTN=TRUE, sample.num=2, SNP.filter=NULL, crv.num=5, cor.threshold=0.3, count.threshold=0.9,
	priority=c("speed", "memory"), bisection.loop=10, ref.gwas=TRUE, theSeed=666666, file.output=TRUE, cpu=10
)
{
#--------------------------------------------------------------------------------------------------------#
# Object: To perform GP/GS with KAML(Genomic Prediction/Selection)										 #
#																										 #
# Input:																								 #
# pfile: phenotype file, one column for a trait, the name of each column must be provided(NA is allowed) #
# gfile: genotype files, including "gfile.geno.desc", "gfile.geno.bin" and "gfile.map"					 #
# kfile: n*n, optional, provided KINSHIP file for all individuals										 #
# dcovfile: n*x, optional, the provided discrete covariates file										 #
# qcovfile: n*x, optional, the provided quantitative covariates file									 #
# pheno：specify phenotype column in the phenotype file(default 1)										 #
# GWAS.model: which model will be used for GWAS(only "GLM" and "MLM" can be selected presently)			 #
# GWAS.npc: the number of PC that will be added as covariance to control population structure			 #
# prior.QTN: the prior QTNs which will be added as covariants, if provided prior QTNs, 	KAML will not	 #
# 			optimize QTNs and model during cross-validation												 #
# prior.model: the prior Model for the prior.QTN that added as covariants								 #
# vc.method: method for variance components estimation("emma" or "gemma")								 #
# K.method: which algorithm will be applied("center", "scale", "vanraden")								 #
# Top.perc: a vector, a subset of top SNPs for each iteration are amplified when calculating KINSHIP	 #
# Top.num: a number, a subset of top SNPs for each iteration are used as covariants						 #
# Logx: a vector, the base for LOG																		 #
# bin.size: the size of each bin																		 #
# max.nQTN: whether limits the max number of Top.num													 #
# sample.num: the sample number of cross-validation														 #
# crv.num: the cross number of cross-validation 														 #
# SNP.filter: the SNPs whose P-value below this threshold will be deleted								 #
# cor.threshold: if the top SNP which is used as covariant is in high correlation with others,			 #
# 			it will be deleted																			 #
# count.threshold: if the count of selected SNP for all iteration >= sample.num*crv.num*count.threshold, #
# 			than it will be treated as covariance in final predicting model								 # 
# bisection.loop: the max loop(iteration) number of bisection algorithm									 #
# ref.gwas: whether to do GWAS for reference population(if not, KAML will merge all GWAS results of	 #
# 			cross-validation by mean)																	 #
# theSeed: the random seed																				 #
# file.output: whether to write the predicted values in file											 #
# cpu: the number of CPU for calculation																 #
#--------------------------------------------------------------------------------------------------------#
	#print the version of KAML
	KAML.version()
    time1 <- as.numeric(Sys.time())
	
	R.ver <- Sys.info()[['sysname']]
	r.open <- !inherits(try(Revo.version, silent=TRUE),"try-error")
	wind <- R.ver == 'Windows'
	linux <- R.ver == 'Linux'
	mac <- (!linux) & (!wind)

	if(wind)	cpu <- 1
	if(r.open &  mac & cpu > 1)	Sys.setenv("VECLIB_MAXIMUM_THREADS" = "1")
	
    #check the parameters
    GWAS.model <- match.arg(GWAS.model)
	K.method <- match.arg(K.method)
	vc.method <- match.arg(vc.method)
	prior.model <- match.arg(prior.model)
	qtn.model <- match.arg(qtn.model)
	priority <- match.arg(priority)
	if(crv.num < 2) stop("'crv.num' must be bigger than 2!")
    if(sample.num < 1) stop("'sample.num' must be bigger than 1!")

	cat(" Attaching data...")
	PHENO <- read.delim(pfile, head=TRUE)
	TAXA <- colnames(PHENO)[pheno]
	PHENO <- PHENO[, pheno]
	N.Ind <- length(PHENO)
	GEBV <- PHENO
	GENO <- attach.big.matrix(paste(gfile, ".geno.desc", sep=""))
	if(length(PHENO) != ncol(GENO))	stop("Number of individuals don't match between pfile and gfile!")
	if(hasNA(GENO@address)){
		stop("NA is not allowed in genotype, use 'KAML.Impute' to impute missings.")
	}
	# MAP <-  try(read.table(paste(gfile, ".map", sep=""), head=FALSE), silent=TRUE)
	# if((!is.null(Top.num) | !is.null(Top.perc)) & class(MAP) == "try-error"){
		# stop("Please provid the Map information for all SNPs!")
	# }
	# if(nrow(MAP) != nrow(GENO)){
		# stop("The number of SNPs in genotype and map doesn't match!")
	# }
	# MAP <- as.matrix(MAP)
	# options(warn = -1)
	# max.chr <- max(as.numeric(MAP[, 2]), na.rm=TRUE)
	# if(is.infinite(max.chr))	max.chr <- 0
	# map.xy.index <- which(!as.numeric(MAP[, 2]) %in% c(0 : max.chr))
	# if(length(map.xy.index) != 0){
		# chr.xy <- unique(MAP[map.xy.index, 2])
		# for(i in 1:length(chr.xy)){
			# MAP[MAP[, 2] == chr.xy[i], 2] <- max.chr + i
		# }
	# }
	# MAP <- matrix(as.numeric(MAP), nrow(MAP))
	# options(warn = 0)
	if(!is.null(kfile)){
		KIN <- read.table(kfile, head=FALSE, colClasses="numeric", stringsAsFactors=FALSE)
		KIN <- as.matrix(KIN)
		if(nrow(KIN) != N.Ind)	stop("Number of individuals don't match between pfile and kfile!")
		if(!is.numeric(KIN)){
			stop("Some none numerical values appear in KINSHIP!")
		}
		if(sum(is.na(KIN))!=0){
			stop("NAs are not allowed in Kinship matrix")
		}
	}
	Cov <- matrix(1, N.Ind)
	if(!is.null(dcovfile)){
		dCov <- read.table(dcovfile, head=FALSE, stringsAsFactors=FALSE)
		if(nrow(dCov) != N.Ind)	stop("Number of individuals don't match between pfile and dcovfile!")
		if(sum(is.na(dCov)) != 0){
			stop("'NA' isn't allowed in 'dcovfile'")
		}
		levelnum <- which(apply(dCov, 2, function(x) length(unique(x))) == nrow(dCov))
		if(length(levelnum) != 0)	stop("No groups in ", paste(levelnum, collapse=","), " column of dcovfile!")
		dCov <- dCov[, apply(dCov, 2, function(x) length(unique(x))) != 1]
		dCov <- as.matrix(dCov)
		X.design <- function(v){
			v=as.character(v)
			vf = factor(v)
			va = as.numeric(vf)
			mrow = length(va)
			mcol = length(levels(vf)) 
			X = matrix(data=c(0),nrow=mrow,ncol=mcol)
			for(i in 1:mrow) {
				ic = va[i]
				X[i,ic] = 1 
			}
			return(X)
		}
		for(i in 1:ncol(dCov)){
			Cov <- cbind(Cov, X.design(dCov[, i])[, -1])
		}
	}
	
	if(!is.null(qcovfile)){
		qCov <- read.table(qcovfile, head=FALSE, stringsAsFactors=FALSE)
		if(nrow(qCov) != N.Ind)	stop("Number of individuals don't match between pfile and qcovfile!")
		if(sum(is.na(qCov)) != 0){
			stop("'NA' isn't allowed in 'qcovfile'")
		}
		Cov <- cbind(Cov, as.matrix(qCov))
	}
	
	cat("Done!\n")
	
	cat(paste(" Trait: ", TAXA, sep=""), "\n")
	NA.Ind <- sum(is.na(PHENO))
	if(NA.Ind == length(PHENO))	stop("No effective values in phenotype file!")
	
	cat(" Number of Total References:", N.Ind-NA.Ind, "\n")
	cat(" Number of Total Predictions:", NA.Ind, "\n")
	cat(" Number of Covariates:", ncol(Cov), "\n")
	cat(" Number of Total SNPs:", nrow(GENO), "\n")
	cat(" Number of CPUs:", cpu, "\n")

	if(is.null(SNP.weight)){
		cat(" New seeds generated from:", theSeed, "\n")
		if((N.Ind-NA.Ind) < 1000 & sample.num < 4 & (!is.null(Top.num) | !is.null(Top.perc))){
			cat(" [Warning: Number of individuals with observations is less than 1000,\n")
			cat("     the predicted GEBV maybe unstable, we recommend setting bigger 'sample.num'!]", "\n")
		}
		if((!is.null(prior.model) && prior.model != "QTN")){
			# cat(" Calculating marker-based Kinship...\n")
			# KIN <- KAML.Kin(GENO, type=K.method, verbose=TRUE); gc()
			if(is.null(kfile)){
				KIN <- KAML.KinNew(GENO, scale=(K.method=="scale"), priority=priority, verbose=TRUE, threads=cpu)
			}
		}else{
			KIN <- NULL
		}

		if(!is.null(prior.QTN)){
			if(is.null(prior.model))	stop("Please provid prior model!")
			Top.num <- NULL
			if(prior.model == "QTN"){Top.num <- NULL; Top.perc <- NULL}
		}else{
			if(prior.model == "QTN")	stop("QTNs must be provided!")
			if(prior.model == "K")	Top.num <- NULL
		}
		cat(" Estimate variance components using:", vc.method,"\n")
		if(!is.null(Top.num) | !is.null(Top.perc)){
			cat(" Performing model: KAML\n")
			cross.res <- KAML.CrossV(phe=PHENO, geno=GENO, prior.QTN=prior.QTN, prior.model=prior.model, vc.method=vc.method, K=KIN, BF.threshold=BF.threshold,
				max.nQTN=max.nQTN, GWAS.model=GWAS.model, theSeed=theSeed, amplify=Logx, K.method=K.method, Top.num=Top.num, binary=binary, priority=priority, 
				Top.perc=Top.perc, sample.num=sample.num, crv.num=crv.num, cpu=cpu, bin.size=bin.size, bisection.loop=bisection.loop, SNP.filter=SNP.filter,
				cor.threshold=cor.threshold, count.threshold=count.threshold, ref.gwas=ref.gwas, GWAS.npc=GWAS.npc, CV=Cov, qtn.model=qtn.model); gc()
			cross.QTN <- cross.res$cross.QTN
			cross.model <- cross.res$cross.model
			cross.tp <- cross.res$cross.tp
			cross.amp <- cross.res$cross.amp
			cross.k <- cross.res$cross.k
			rm(cross.res); gc()
		}else{
			if(is.null(prior.QTN)){
				cat(" Performing model: GBLUP\n")
				cross.QTN <- NULL
				cross.model <- "K"
				cross.tp <- NULL
				cross.amp <- NULL
				cross.k <- KIN
			}else{
				cat(" Performing model: KAML\n")
				cat(" The provided QTNs:", prior.QTN, "\n")
				cat(paste(" The provided Model: ", prior.model, "\n", sep=""))
				cross.QTN <- prior.QTN
				cross.model <- prior.model
				cross.tp <- NULL
				cross.amp <- NULL
				cross.k <- KIN
			}
		}
		rm("KIN"); gc()
		cat(" Predicting...\n")
		if(is.null(cross.QTN)){
			myest <- KAML.Mix(phe=PHENO, K=cross.k, vc.method=vc.method, CV=Cov)
			GEBV <- myest$ebv
			beta <- as.vector(myest$beta)
		}else
		{
			if(cross.model == "QTN+K"){
				QTN.cv <- cbind(Cov, t(GENO[cross.QTN, , drop=FALSE]))
				myest <- KAML.Mix(phe=PHENO, CV=QTN.cv, vc.method=vc.method, K=cross.k)
				GEBV <- myest$ebv + t(GENO[cross.QTN, , drop=FALSE]) %*% as.matrix(myest$beta[-c(1:ncol(Cov))])
				beta <- as.vector(myest$beta[c(1:ncol(Cov))])
			}
			if(cross.model == "QTN"){
				QTN.cv <- t(GENO[cross.QTN, , drop=FALSE])
				NA.index <- is.na(PHENO)
				qtn.eff <- KAML.GLM(y=PHENO[!NA.index], X=Cov[!NA.index, ], qtn.matrix=QTN.cv[!NA.index, ])$QTN.eff
				GEBV <- data.matrix(QTN.cv) %*% as.matrix(qtn.eff[-c(1:ncol(Cov))])
				beta <- as.vector(qtn.eff[c(1:ncol(Cov))])
			}
		}
	}else{
		if(length(SNP.weight) != nrow(GENO))	stop("length of weights should equal to the number of SNPs in genotype!")
		# K <- KAML.Kin(GENO, type="center", priority, weight = SNP.weight, SUM=nrow(geno))
		K <- KAML.KinNew(GENO, SUM=nrow(geno), scale=(K.method=="scale"), priority=priority, weight = SNP.weight, verbose=TRUE, threads=cpu)
			
		myest <- KAML.Mix(phe=PHENO, CV=Cov, vc.method=vc.method, K=K)
		GEBV <- myest$ebv
		beta <- as.vector(myest$beta[c(1:ncol(Cov))])
		cross.QTN <- NULL
		cross.model <- "K"
		cross.tp <- NULL
		cross.amp <- NULL
		cross.k <- K
	}
	#return the results
	GEBV <- as.matrix(GEBV)
	colnames(GEBV) <- TAXA
	yHat <- Cov %*% as.matrix(beta) + GEBV

	if(file.output == TRUE){
		file.name <- paste("KAML.", TAXA, ".pred", ".txt", sep="")
		write.table(yHat, file.name, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
	}
	
	time2 <- as.numeric(Sys.time())
	time.cal <- round(time2-time1)
	times <- function(x){
		h <- x %/% 3600
		m <- (x %% 3600) %/% 60
		s <- ((x %% 3600) %% 60)
		index <- which(c(h, m, s) != 0)
		num <- c(h, m, s)[index]
		char <- c("h", "m", "s")[index]
		return(paste(round(num), char, sep="", collapse=""))
	}
	cat(paste(" ", TAXA, " is DONE within total run time: ", times(time.cal), "\n", sep=""))
	cat(paste(c("#", rep("-", 19), "KAML ACCOMPLISHED SUCCESSFULLY", rep("-", 19), "#"), collapse=""),"\n\r")
	return(list(y=yHat, beta=beta, gebv=GEBV, qtn=cross.QTN, model=cross.model, top.perc=cross.tp, logx=cross.amp, K=cross.k))
}

KAML.EIGEN.REML <- 
function(y, X, eigenK)
{
	reml <- lmm.diago(Y=y, X=X, eigenK=eigenK, method="brent", verbose=FALSE)
	vg <- reml[[2]]
	ve <- reml[[1]]
	delta <- ve / vg
	return(list(vg=vg, ve=ve, delta=delta))
}

Math_cpu_check <- 
function()
{
    r.open <- !inherits(try(Revo.version, silent=TRUE),"try-error")
	if(r.open){
		cpu <- try(getMKLthreads(), silent=TRUE)
		if(class(cpu) == "try-error"){
			return(2)
		}else{
			return(cpu)
		}
	}else{
		return(1)
	}
}
