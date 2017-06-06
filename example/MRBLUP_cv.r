MRBLUP.version <-
function()
{
#--------------------------------------------------------------------------------------------------------#
# KAMLMM: Kinship Adjusted Multiple-locus Linear Mixed Model											 #
# Writen by Lilin Yin and Xiaolei Liu																	 #
# Last update: May 20, 2017																				 #
#--------------------------------------------------------------------------------------------------------#
MRBLUP.versions <- "1.0.1"
cat(paste("#", paste(rep("-", 27), collapse=""), "Welcome to KAMLMM", paste(rep("-", 27), collapse=""), "#", sep=""), "\n")
cat("#", paste(rep(" ", 69), collapse=""), "#", "\n")
cat("#", rep(" ", 5), "Kinship Adjusted Multiple-locus Linear Mixed Model ", rep(" ", 4), "#", "\n")
cat("#", paste(rep(" ", 69), collapse=""), "#", "\n")
cat(paste("#  ", "Version: ", MRBLUP.versions, paste(rep(" ", 55), collapse=""), "#", sep=""), "\n")
cat(paste("#  ", "Authors: Lilin Yin & Xiaolei Liu", paste(rep(" ", 24), collapse=""), "_\\\\|//_", paste(rep(" ", 6), collapse=""), "#", sep=""), "\n")
cat(paste("#  ", "Contact: ylilin@163.com & xll19870827@hotmail.com", paste(rep(" ", 6), collapse=""), "//^. .^\\\\",
	paste(rep(" ", 5), collapse=""), "#", sep=""), "\n")
cat(paste("#", paste(rep("-", 53), collapse=""), "ooO-( (00) )-Ooo", paste(rep("-", 2), collapse=""), "#", sep=""), "\n")
}

BayesR <- 
function(
	bfile="", n=1, seed=123456, file.remove=TRUE, memo='BayesR'
)
{
	cat("------------BayesR started------------\n")
	phe <- read.delim(paste(bfile, ".fam", sep=""), head=FALSE, sep=" ")
	map <- read.delim(paste(bfile, ".bim", sep=""), head=FALSE, sep="\t")
	indexNA <- which(is.na(phe[, 5+n]))
	cat(paste("Number of Individuals:", nrow(phe), "\n"))
	cat(paste("Number of References:", nrow(phe)-length(indexNA), "\n"))
	cat(paste("Number of Candidates:", length(indexNA), "\n"))
	cat(paste("Number of SNPs:", nrow(map), "\n"))
	rm(phe);rm(map);gc()
	
	cat("Running...\n")
	
	t1=as.numeric(Sys.time())
	
	system(
		paste(
			"bayesR -bfile", bfile,
			"-n", n,
			"-out", memo,
			"-seed", seed
		)
	)
	
	system(
		paste(
			"bayesR -bfile", bfile,
			"-n", n,
			"-out", paste(memo, "pred", sep="_"),
			"-predict -model", paste(memo, ".model", sep=""),
			"-freq", paste(memo, ".frq", sep=""),
			"-param", paste(memo, ".param", sep="")
		)
	)
	
	t2=as.numeric(Sys.time())
	time <- t2 - t1
	
	refebv <- as.matrix(read.table(paste(memo, ".gv", sep=""), head=FALSE)[, 1])
	infebv <- as.matrix(read.delim(paste(paste(memo, "pred", sep="_"), ".gv", sep=""), head=FALSE)[, 1])
	if(length(indexNA) == 0){
		refebv <- refebv
		infebv <- NULL
	}else{
		refebv <- refebv[-indexNA, ]
		infebv <- infebv[indexNA, ]
	}
	if(file.remove){
		files <- paste(memo, c(".model", ".frq", ".param", ".hyp", ".gv", ".log", "_pred.gv", "_pred.log"), sep="")
		unlink(files)
	}
	cat("-----------------DONE-----------------\n")
	return(list(refebv=refebv, infebv=infebv, time=time))
}

BSLMM <- 
function(
	bfile="", n=1, seed=123456, file.remove=TRUE, memo='BSLMM'
)
{
	cat("-------------BSLMM started------------\n")
	phe <- read.delim(paste(bfile, ".fam", sep=""), head=FALSE, sep=" ")
	map <- read.delim(paste(bfile, ".bim", sep=""), head=FALSE, sep="\t")
	indexNA <- which(is.na(phe[, 5+n]))
	cat(paste("Number of Individuals:", nrow(phe), "\n"))
	cat(paste("Number of References:", nrow(phe)-length(indexNA), "\n"))
	cat(paste("Number of Candidates:", length(indexNA), "\n"))
	cat(paste("Number of SNPs:", nrow(map), "\n"))
	rm(phe);rm(map);gc()
	
	cat("Running...\n")
	
	t1=as.numeric(Sys.time())
	
	system(
		paste(
			"gemma -bfile", bfile,
			"-n", n,
			"-notsnp -bslmm 1 ",
			"-o", memo,
			"-seed", seed,
			">", paste(memo, ".LOG", sep="")
		)
	)
	
	File <- read.delim(paste("./output/", memo, ".log.txt", sep=""), head=FALSE)
	cat(paste("Number of analyzed SNPs:", as.numeric(unlist(strsplit(as.character(File[14,]), " "))[7]), "\n"))
	
	if(length(indexNA) != 0){
		system(
			paste(
				"gemma -bfile ", bfile,
				" -n ", n,
				" -epm ./output/", memo, ".param.txt",
				" -emu ./output/", memo, ".log.txt",
				" -seed ", seed,
				" -predict 1 -o ",paste(memo ,"pred",sep="_"),
				" >>", paste(memo, ".LOG", sep="")
			,sep="")
		)
	}
	t2=as.numeric(Sys.time())
	time <- t2 - t1
	
	refebv <- as.matrix(read.table(paste("./output/", memo, ".bv.txt", sep=""), head=FALSE)[, 1])

	if(length(indexNA) == 0){
		refebv <- refebv
		infebv <- NULL
	}else{
		infebv <- as.matrix(read.table(paste("./output/", paste(memo, "pred", sep="_"), ".prdt.txt", sep=""), head=FALSE)[, 1])
		refebv <- refebv[-indexNA, ]
		infebv <- infebv[indexNA, ]
	}
	if(file.remove){
		unlink(
			c(paste(memo, ".LOG", sep=""),
			paste("./output/",memo,
			c(".bv.txt",".gamma.txt",".hyp.txt",".log.txt",".param.txt",
			"_pred.log.txt","_pred.prdt.txt"),sep="")
			)
		)
	}
	cat("-----------------DONE-----------------\n")
	return(list(refebv=refebv, infebv=infebv, time=time))
}

MultiBLUP <- 
function(
	bfile='', n=1, wind=75000, sig1=NULL, sig2=0.01, background=TRUE, file.remove=TRUE, memo='multiBLUP'
)
{

	cat("----------------------MultiBLUP started----------------------\n")
	phe <- read.delim(paste(bfile, ".fam", sep=""), head=FALSE, sep=" ")
	map <- read.delim(paste(bfile, ".bim", sep=""), head=FALSE, sep="\t")
	indexNA <- which(is.na(phe[, 5+n]))
	phe <- phe[!is.na(phe[, 5+n]), ]
	cat(paste("Number of Individuals:", nrow(phe), "\n"))
	cat(paste("Number of References:", nrow(phe)-length(indexNA), "\n"))
	cat(paste("Number of Candidates:", length(indexNA), "\n"))
	cat(paste("Number of SNPs:", nrow(map), "\n"))

	dir.create(memo)
	write.table(phe[,c(1,2,5+n)], paste("./", memo, "/", memo, ".pheno", sep=""), row.names=FALSE, col.names=FALSE, quote=FALSE, sep=" ")
	rm(phe);rm(map);gc()

	t1 <- as.numeric(Sys.time())
	
	#divide the genome into chunks. 
	cat(paste("Divide the genome into chunks(", wind/1000, "kb)...", sep=""))
	multiBLUP.try <- try(
		system(
			paste(
				"ldak --cut-genes ", memo, 
				" --chunks-bp ", wind, 
				" --bfile ", bfile, 
				" --ignore-weights YES",
				" >", paste(memo, "/", memo, "_multiBLUP.log", sep=""),
			sep="")
		), silent=TRUE
	)
	cat("Done!\n")

	#test the chunks for association.
	cat("Test the chunks for association...")
	multiBLUP.try <- try(
		system(
				paste(
				"ldak --calc-genes-reml ", memo,
				" --bfile ", bfile,
				" --minmaf 0.01 --minobs 0.90 " ,
				" --pheno ", paste(memo, "/", memo, ".pheno", sep=""),
				" --ignore-weights YES --partition 1",
				" >>", paste(memo, "/", memo, "_multiBLUP.log", sep=""),
			sep="")
		), silent=TRUE
	)
	cat("Done!\n")

	#identify significant chunks.
	if(is.null(sig1)){
		chunks.details <- read.delim(paste(memo, "/regress1", sep=""), head=TRUE, sep=" ")
		sig1 <- 0.05 / nrow(chunks.details)
		rm(chunks.details); gc()
	}
	cat(paste("Identify significant chunks(", sig1, " ", sig2, ")...", sep=""))
	multiBLUP.try <- try(
		system(
			paste(
				"ldak --join-genes-reml ", memo,
				" --bfile ", bfile,
				" --sig1 " , sig1,
				" --sig2 ", sig2,
				" >>", paste(memo, "/", memo, "_multiBLUP.log", sep=""),
			sep="")
		), silent=TRUE
	)
	cat("Done!\n")

	regionN <- read.delim(paste("./", memo, "/region_number", sep=""), head=FALSE)
	cat(paste("Number of regions:", as.numeric(regionN), "\n"))
	
	if(regionN == 0 | background){
		#calculate the background kinship matrix.
		cat("Calculate the background kinship matrix...")
		multiBLUP.try <- try(
			system(
				paste(
					"ldak --calc-kins-direct ", paste(memo, "/region0", sep=""),
					" --bfile ", bfile,
					" --extract " , paste(memo, "/region0", sep=""),
					" --minmaf 0.01 --minobs 0.90 "  ,
					" --ignore-weights YES",
					" >>", paste(memo, "/", memo, "_multiBLUP.log", sep=""),
				sep="")
			), silent=TRUE
		)
		cat("Done!\n")
	}
	
	#estimate variance components for the MultiBLUP Model.
	cat("Estimate variance components for the MultiBLUP Model...")
	if(as.numeric(regionN) != 0 & !background){
		multiBLUP.try <- try(
			system(
				paste(
					"ldak --reml ", paste(memo, "/REML", sep=""),
					" --bfile ", bfile,
					" --pheno ", paste(memo, "/", memo, ".pheno", sep=""),
					" --region-number " , as.numeric(regionN),
					" --region-prefix " , paste(memo, "/region", sep=""),
					" --ignore-weights YES",
					" >>", paste(memo, "/", memo, "_multiBLUP.log", sep=""),
				sep="")
			), silent=TRUE
		)
	}else{
		multiBLUP.try <- try(
			system(
				paste(
					"ldak --reml ", paste(memo, "/REML", sep=""),
					" --bfile ", bfile,
					" --grm " , paste(memo, "/region0", sep=""),
					" --pheno ", paste(memo, "/", memo, ".pheno", sep=""),
					" --region-number " , as.numeric(regionN),
					" --region-prefix " , paste(memo, "/region", sep=""),
					" --ignore-weights YES",
					" >>", paste(memo, "/", memo, "_multiBLUP.log", sep=""),
				sep="")
			), silent=TRUE
		)
	}
	cat("Done!\n")

	#estimate SNP effect sizes and calculate predictions.
	cat("Estimate SNP effect sizes and calculate predictions...")
	if(as.numeric(regionN) != 0 & !background){
		multiBLUP.try <- try(
			system(
				paste(
					"ldak --calc-blups ", paste(memo, "/PRED", sep=""),
					" --bfile ", bfile,
					" --remlfile ", paste(memo, "/REML.reml", sep=""),
					" >>", paste(memo, "/", memo, "_multiBLUP.log", sep=""),
				sep="")
			), silent=TRUE
		)
	}else{
		multiBLUP.try <- try(
			system(
				paste(
					"ldak --calc-blups ", paste(memo, "/PRED", sep=""),
					" --bfile ", bfile,
					" --grm " , paste(memo, "/region0", sep=""),
					" --remlfile ", paste(memo, "/REML.reml", sep=""),
					" >>", paste(memo, "/", memo, "_multiBLUP.log", sep=""),
				sep="")
			), silent=TRUE
		)
	}
	cat("Done!\n")
	
	t2 <- as.numeric(Sys.time())
	time <- t2 - t1
	
	cat("----------------------------DONE-----------------------------\n")
	if(class(multiBLUP.try) == "try-error"){
		return("error!")
	}else{
		gebv <- as.matrix(read.delim(paste(memo,"/PRED.pred",sep=""), head=TRUE)[, 3])
		if(length(indexNA) == 0){
			refebv <- gebv
			infebv <- NULL
		}else{
			refebv <- gebv[-indexNA, ]
			infebv <- gebv[indexNA, ]
		}
		if(file.remove){
			unlink(memo, recursive=TRUE)
		}
		return(list(refebv=refebv, infebv=infebv, time=time))
	}
}

BayesABCL <-
function(
	refphe, refgeno, infgeno, Model=c("BayesA", "BayesB", "BayesC", "BayesLASSO"), memo="", niter=100000, burnIn=20000, thin=1
)
{
	memo <- paste(memo, Model, ".", sep="")
	packages <- rownames(installed.packages())
	if("BGLR" %in% packages){
		print("BGLR is available!")
	}else
	{
		print("BGLR is not available!")
		print("Installing the BGLR package from CRAN!")
		install.packages(pkg='BGLR', repos='https://cran.r-project.org/')
	}
	library("BGLR")
    geno <- cbind(refgeno, infgeno)
	#print("Normalizing the genotype...")
	geno <- t(geno);gc()
    #geno <- scale(geno)
    #geno <- geno[ ,!is.na(geno[1,])]
	#print("Normalizing DONE!")
    n.inf <- ncol(infgeno)
	y <- c(as.numeric(as.character(refphe[, 2])), rep(NA, n.inf))
	rm(list=c("refgeno", "infgeno"))
	gc()
    time1 <- as.numeric(Sys.time())
	print("BGLR started...")
	switch(
		match.arg(Model), 
		"BayesA"={
			ETA<-list(list(X=geno, model='BayesA'))
			blm=BGLR(y=y, ETA=ETA, nIter=niter, burnIn=burnIn, thin=thin, saveAt=memo, df0=5, S0=NULL)
		}, 
		"BayesB"={
			ETA<-list(list(X=geno, model='BayesB', probIn=0.05))
			blm=BGLR(y=y, ETA=ETA, nIter=niter, burnIn=burnIn, thin=thin, saveAt=memo, df0=0.00429, S0=4.234)
		}, 
		"BayesC"={
			ETA<-list(list(X=geno, model='BayesC'))
			blm=BGLR(y=y, ETA=ETA, nIter=niter, burnIn=burnIn, thin=thin, saveAt=memo, df0=5, S0=NULL)
		}, 
		"BayesLASSO"={
			ETA<-list(list(X=geno, model='BL'))
			blm=BGLR(y=y, ETA=ETA, nIter=niter, burnIn=burnIn, thin=thin, saveAt=memo, df0=5, S0=NULL)
		}
	)
    time2 <- as.numeric(Sys.time())
    timer <- time2-time1
	rm("geno")
	gc()
	if(Model == "BayesLASSO"){
		unlink(paste(memo,c("varE.dat","mu.dat","ETA_1_lambda.dat"),sep=""))
	}else{
		unlink(paste(memo,c("varE.dat","mu.dat",paste("ETA_1_par",Model,".dat",sep="")),sep=""))
	}
	snp.eff <- blm$ETA[[1]]$b
	gebv <- tail(blm$yHat, n=n.inf)
	return(list(snp.eff=snp.eff, gebv=gebv, timer=timer))
}

MRBLUP.stas.cal <- 
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

MRBLUP.Bar <- 
function(
	i, n, type=c("type1", "type2", "type3"), symbol="-", tmp.file=NULL, symbol.head="|", symbol.tail=">" ,fixed.points=TRUE, points=seq(0,100,5), symbol.len=50
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
							sprintf("%.2f%%", 100*i/n), "\n", sep="")
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
						sprintf("%.2f%%", 100*i/n), "\n", sep="")
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

MRBLUP.File <- 
function(
	phe, geno, map, memo="mydata", type=c("BIMBAM", "Binary","PED")
)
{
#--------------------------------------------------------------------------------------------------------#
# Object: transform numeric genotype into "BIMBAM", "Bianry", "PED" 									 #
#	 																									 #
# Input:	 																							 #
# phe: a numeric vector(NA is allowed)	 																 #
# geno: m*n matrix(0, 1, 2), m is the number of markers, n is the number of individuals					 #
# map: m*3 matrix, the three columns are SNP, Chr, Pos respectively	 									 #
# memo: the name of output files	 																	 #
# type: which file type will be transformed into	 													 #
#--------------------------------------------------------------------------------------------------------#
	phe[is.na(phe)]=0
	geno=as.matrix(geno)
	map=as.matrix(map)
	map=cbind(as.character(map[,1]),as.integer(map[,2]),as.integer(map[,3]))
	switch(
	match.arg(type),
	"PED"={
		cat("Transpose...")
		geno=t(geno); gc()
		cat("Done!\n")
		cat("Replace...")
		geno[is.na(geno)]=00
		geno[geno == 0]=11
		geno[geno == 1]=12
		geno[geno == 2]=22
		cat("Done!\n")
		file.map <- cbind(map[, 2], map[, 1], 0, map[, 3])
		file.ped <- cbind(c(1:length(phe)), c(1:length(phe)), matrix(0, length(phe), 3), phe, geno)
		cat("Writing map...")
		write.table(file.map, paste(memo, ".map", sep=""), col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
		cat("Done!\n")
		cat("Writing ped...\n")
		fileNumCon<-file(description=paste(memo,".ped",sep=""), open="w")
		for(i in 1:nrow(file.ped)){
			MRBLUP.Bar(i = i, n = nrow(file.ped))
			writeLines(as.character(file.ped[i,1:(ncol(file.ped)-1)]), fileNumCon, sep="\t")
			writeLines(as.character(file.ped[i,ncol(file.ped)]), fileNumCon, sep="\n")
		}
		close.connection(fileNumCon)
		cat("DONE!\n")
	},
	"BIMBAM"={
		pheno=as.matrix(phe)
		geno=cbind(map[, 1], 0, 0, geno)
		anno=cbind(map[, 1], map[, 3], map[, 2], 0)
		cat("Writing phenotype...")
		write.table(pheno, paste(memo, ".pheno.txt", sep=""), col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
		cat("Done!\n")
		cat("Writing genotype...")
		write.table(geno, paste(memo, ".geno.txt", sep=""), col.names=FALSE, row.names=FALSE, quote=FALSE, sep=", ")
		cat("Done!\n")
		cat("Writing map...")
		write.table(anno, paste(memo, ".anno.txt", sep=""), col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
		cat("Done!\n")
	}, 
	"Binary"={
		files <- list.files()
		cat("Transpose...")
		geno=t(geno)
		cat("Done!\n")
		cat("Replace...")
		geno[geno == 0]=11
		geno[geno == 1]=12
		geno[geno == 2]=22
		cat("Done!\n")
		file.map <- cbind(map[, 2], map[, 1], 0, map[, 3])
		file.ped <- cbind(c(1:length(phe)), c(1:length(phe)), matrix(0, length(phe), 3), phe, geno)
		cat("Writing map...")
		write.table(file.map, paste(memo, ".map", sep=""), col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
		cat("Done!\n")
		cat("Writing ped...\n")
		fileNumCon<-file(description=paste(memo,".ped",sep=""), open="w")
		for(i in 1:nrow(file.ped)){
			MRBLUP.Bar(i = i, n = nrow(file.ped))
			writeLines(as.character(file.ped[i,1:(ncol(file.ped)-1)]), fileNumCon, sep="\t")
			writeLines(as.character(file.ped[i,ncol(file.ped)]), fileNumCon, sep="\n")
		}
		close.connection(fileNumCon)
		cat("Plinking...\n")
		system(paste("plink --file ", memo,
			" --allow-extra-chr --chr-set ", length(unique(map[,2]))-4,
			" --compound-genotypes --out ", memo,
			" --make-bed", " > /dev/null 2>&1", sep=""))
		cat("DONE!\n")
		filesN <- list.files()
		unlink(filesN[!filesN %in% c(files, paste(memo, c(".bed", ".bim", ".fam", ".log"), sep=""))], recursive = TRUE)
	}
	)
}

MRBLUP.Bin <-
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

MRBLUP.GEMMA.VC <-
function(
	y, X=NULL, K, memo="trait"
)
{
#--------------------------------------------------------------------------------------------------------#
# Object: To estimate variance component using GEMMA software											 #
#	 																									 #
# Input:	 																							 #
# y: The value vector	 																				 #
# X: The fixed effect(X must contain a column of 1's)													 #
# K: Kinship for all individuals(row and column orders must be the same with y)	 						 #
# memo: A character that will be added to the names of output files										 #
#--------------------------------------------------------------------------------------------------------#
    if(system("gemma > /dev/null 2>&1") == 127)   stop("Please download GEMMA and configure it to the environment path")

	y=as.matrix(y)
	if(!is.null(X) & (max(X)!=min(X))){
		write.table(X, paste(memo, ".cv.txt", sep=""), col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
	}
	write.table(y, paste(memo, ".pheno.txt", sep=""), col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
	write.table(K, paste(memo, ".Kin.txt", sep=""), col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
	if(!is.null(X) & (max(X)!=min(X))){
		system(paste("gemma -p", paste(memo, ".pheno.txt", sep=""), "-c", paste(memo, ".cv.txt", sep=""), "-k",
        paste(memo, ".Kin.txt", sep=""), "-n 1 -vc 2 -o", memo, "> /dev/null 2>&1"))
	}else
	{
		system(paste("gemma -p", paste(memo, ".pheno.txt", sep=""), "-k", paste(memo, ".Kin.txt", sep=""),
        "-n 1 -vc 2 -o", memo, "> /dev/null 2>&1"))
	}
	gemma.res=read.delim(paste("./output/", memo, ".log.txt", sep=""), head=FALSE)
	vg=as.numeric(unlist(strsplit(as.character(gemma.res[18, ]), "  "))[2])
	ve=as.numeric(unlist(strsplit(as.character(gemma.res[18, ]), "  "))[3])
	delta=round(ve/vg, 4)
	unlink(c(paste(memo, ".pheno.txt", sep=""), paste(memo, ".Kin.txt", sep=""), paste("./output/", memo, ".log.txt", sep="")), recursive = TRUE)
	if(!is.null(X)) unlink(paste(memo, ".cv.txt", sep=""), recursive = TRUE)
	return (list(delta=delta, ve=ve, vg=vg))
}

MRBLUP.GEMMA.REML <- 
function(
	y, X, K, rtol = 1e-6, atol = 1e-8, ctol = 1e-8, root = FALSE
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
# atol and ctol: the desired accuracy																	 #
# root: whether to optimize parameters by multiroot(the will be slow)			 						 #
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
		WtW = crossprod(X)
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
		d = crossprod(Kry)
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
			para2 = (-0.5 * sum(diag(P)) + 0.5 * crossprod(Py)) * exp(log_sigma2[2]))
	}

    # multiroot
	if(root){
		mult.res = rootSolve::multiroot(LogRL_dev1, start = log_sigma2, parms=list(y = y,K = K,X = X),rtol = rtol, atol = atol, ctol = ctol, maxiter=1000)
		vg = exp(mult.res$root[1])
		ve = exp(mult.res$root[2])
		delta = ve / vg
	}else{
		vg = exp(log_sigma2[1])
		ve = exp(log_sigma2[2])
		delta = ve / vg
	}
	return(list(vg=vg, ve=ve, delta=delta))
}

MRBLUP.Mix <- 
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
	inf.index <- which(is.na(phe))
	ref.index <- c(1:length(phe))[-inf.index]
	refphe <- phe[ref.index]
	N <- length(phe)
	n <- length(refphe)
	
	if(is.null(K)){
		K <- MRBLUP.Kin(geno)
	}
	
	if(is.null(CV)){
		X <- matrix(1, n, 1)
	}else{
		X <- as.matrix(CV[ref.index, ])
	}

	if(is.null(lambda)){
        if(vc.method == "emma") {
			emma.reml <- MRBLUP.EMMA.REML(refphe, X=X, K=K[ref.index, ref.index])
			lambda <- emma.reml$delta
			emma.ll <- emma.reml$REML
		}
		if(vc.method == "gemmaU") lambda <- MRBLUP.GEMMA.REML(refphe, X=X, K=K[ref.index, ref.index], root=FALSE)$delta
	}

	#there is a error when running in Mcrosoft R open with parallel
	if(!is.null(eigen.K)){
		eig <- eigen.K
	}else{
		math.cpu <- try(getMKLthreads(), silent=TRUE)
		if(length(ref.index) < 2000) {
            try(setMKLthreads(1), silent=TRUE)
        }
        eig <- eigen((K[ref.index, ref.index]), symmetric=TRUE)
        if(length(ref.index) < 2000) {
            try(setMKLthreads(math.cpu), silent=TRUE)
        }

	}
	y <- as.numeric(as.character(refphe))
	U <- eig$vectors * matrix(sqrt(1/(eig$values + lambda)), n, length(eig$values), byrow=TRUE)
	yt <- crossprod(U, y)
    X0t <- crossprod(U, X)
	Xt <- X0t
	X0X0 <- crossprod(X0t)
	if(X0X0[1, 1] == "NaN"){
		Xt[which(Xt == "NaN")]=0
		yt[which(yt == "NaN")]=0
		XX=crossprod(Xt)
	}
	X0Y <- crossprod(X0t, yt)
	XY <- X0Y
	iX0X0 <- try(solve(X0X0), silent = TRUE)
	if(inherits(iX0X0, "try-error")){
		if(r.open){
			iX0X0 <- MASS::ginv(X0X0)
		}else{
			iX0X0 <- try(geninv(X0X0), silent = TRUE)
			if(inherits(iX0X0, "try-error")){
				warning("Package 'rfunctions' is not installed!")
				iX0X0 <- MASS::ginv(X0X0)
			}
		}	
	}
	iXX <- iX0X0
	beta <- crossprod(iXX, XY)
	XtimesBetaHat <- X %*% beta
	YminusXtimesBetaHat <- matrix(y)- XtimesBetaHat
	Dt <- crossprod(U, YminusXtimesBetaHat)
	Zt <- t(U); rm(U); gc()
	if(X0X0[1, 1] == "NaN"){
		Dt[which(Dt == "NaN")]=0
		Zt[which(Zt == "NaN")]=0
	}
	Z <- diag(N)
	Z[inf.index, inf.index] <- 0
	Z <- Z[ref.index, ]
    KZt <- tcrossprod(K, Z)
    S <- crossprod(Zt, Dt)
	BLUP.ebv <- as.vector(KZt %*% S)
	return(list(beta=beta, ebv=BLUP.ebv, LL=emma.ll))
}

MRBLUP.Kin <- 
function(
	M, weight=NULL, type="center", priority=c("speed", "memory"), memo=NULL, SUM=NULL, maxLine=1000
)
{
#--------------------------------------------------------------------------------------------------------#
# Object: To calculate the Vanraden Kinship																 #
#																										 #
# Input:	 																							 #
# M: Genotype in numeric format, pure 0, 1, 2 matrix; m * n, m is marker size, n is population size	 	 #
# (both matrix or big.matrix are allowed for M)															 #
# weight: the weights for all makers, the length of weight is m											 #
# type: which algorithm will be applied("center", "scale", "vanraden")									 #
# priority: choose the calculation speed or memory														 #
# memo: the names of temporary files																	 #
# SUM: the sum will be used to scale the matrix						 									 #
# maxLine：this parameter can control the memory size when using big.matrix								 #
#--------------------------------------------------------------------------------------------------------#
	if(!type %in% c("scale", "center", "vanraden"))	stop("please select the right kinship algorithm: 'center', 'scale', 'vanraden'!")
	if(!is.null(weight)){
		if(sum(is.na(weight)) != 0)	stop("'NA' is not allowed in weight")
	}
    if(is.null(dim(M))) M <- t(as.matrix(M))
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
			if(!check.r){
				#rfunctions package must be installed
				K <- try(crossprodcpp(M)/SUM, silent=TRUE)
				if(inherits(K, "try-error")){
					warning("Package 'rfunctions' is not installed!")
					K <- crossprod(M)/SUM
				}
			}else{
				#crossprod using parallel in microsoft r open
				K <- 0.5 * crossprod(M)/SUM
			}	
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
			K <- big.crossprod(Z)/SUM
			rm(Z); gc()
			unlink(c(paste("Z", memo, ".temp.bin", sep=""), paste("Z", memo, ".temp.desc", sep="")), recursive = TRUE)
		}
	)
	return(K)
}

MRBLUP.EMMA.REML <- 
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
		cXX <- try(crossprodcpp(X), silent = TRUE)
		if(inherits(cXX, "try-error")){
			warning("Package 'rfunctions' is not installed!")
			cXX <- try(crossprod(X), silent = TRUE)
		}
		iXX <- try(solve(cXX), silent = TRUE)
		if(inherits(iXX, "try-error")){
			if(r.open){
				iXX <- MASS::ginv(cXX)
			}else{
				iXX <- try(geninv(cXX), silent = TRUE)
				if(inherits(iXX, "try-error")){
					warning("Package 'rfunctions' is not installed!")
					iXX <- MASS::ginv(cXX)
				}
			}
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
		if(ncol(K) < 2000) {
			try(setMKLthreads(1), silent=TRUE)
		}
        eig <- eigen(S %*% (K + diag(1, n)) %*% S, symmetric=TRUE)#S4 error here
		if(ncol(K) < 2000) {
			try(setMKLthreads(math.cpu), silent=TRUE)
		}
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
	cXX <- try(crossprodcpp(X), silent=TRUE)
	if(inherits(cXX, "try-error")){
		warning("Package 'rfunctions' is not installed!")
		cXX <- crossprod(X)
	}
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

MRBLUP.Data <- 
function(
	hfile=NULL, numfile=NULL, mapfile=NULL, bfile=NULL, out=NULL, sep="\t", SNP.impute="Middle", maxLine=10000, priority="memory"
)
{
#--------------------------------------------------------------------------------------------------------#
# Object: To prepare data for MRBLUP package															 #
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
	cat("Preparing data for MRBLUP...\n")
	#if(is.null(hfile)&is.null(numfile))
	#stop("Hapmap or Numeric genotype is needed.")
	if(!is.null(hfile)&!is.null(numfile))
	stop("Only one of Hapmap or Numeric genotype is needed!")
	if((!is.null(numfile) & is.null(mapfile)) | (is.null(numfile) & !is.null(mapfile)))
	stop("Both Map and Numeric genotype are needed!")
	if(is.null(out)) out="MRBLUP"

	if(!is.null(bfile)){
		map <- read.table(paste(bfile, ".bim", sep=""), head=FALSE)
		if(length(unique(map[, 2])) != nrow(map)){
			stop("The names of all SNPs must be unique, please check 'BIM' file!")
		}
		map <- map[, c(2, 1, 4)]
		cat("Reading binary files...")
		write.table(map, paste(out, ".map", sep=""), row.names=FALSE, col.names=FALSE, sep=sep, quote=FALSE)
		geno <- read.plink(paste(bfile, ".bed", sep=""))[[1]]
		cat("Done!\n")
		bck <- paste(out, ".geno.bin", sep="")
		dsc <- paste(out, ".geno.desc", sep="")
		nmarkers <- ncol(geno)
		ns <- nrow(geno)
		myGeno.backed<-big.matrix(nmarkers, ns, type="char",
			backingfile=bck, descriptorfile=dsc)
		options(bigmemory.typecast.warning=FALSE)
		cat("Output BIG genotype...\n")
		inGENOFile=TRUE
		i <- 0
		#printN <- unlist(strsplit(x=as.character(nmarkers), split="", fixed=TRUE))
		#printIndex <- seq(0, (as.numeric(printN[1]) + 1) * 10^(length(printN)), 1000)[-1]
		Num.fun <- function(x){
			x <- data.matrix(as.data.frame(x))
			x[x==0]=2
			return(x)
		}
		while(inGENOFile){
			i <- i + maxLine
			if(i >= nmarkers){
				xx <- nmarkers
				inGENOFile <- FALSE
			}else{
				xx <- i
			}
			#if(sum(i >= printIndex )!=0){
			#	printIndex <- printIndex[printIndex > i]
			#	print(paste("Number of Markers Written into BIG File: ", xx, sep=""))
			#}
			if(i >= nmarkers){
				myGeno.backed [(i-maxLine + 1):nmarkers, ] <- -1 * apply(geno[, (i-maxLine + 1):nmarkers], 1, Num.fun) + 3
			}else{
				myGeno.backed [(i-maxLine + 1):i, ] <- -1 * apply(geno[, (i-maxLine + 1):i], 1, Num.fun) + 3
			}
			MRBLUP.Bar(i=xx, n=nmarkers, fixed.points=FALSE)
		}
		geno.flush <- flush(myGeno.backed)
		if(!geno.flush){
			stop("flush failed")
		}else{
			cat("Preparation for BIG data is done!\n")
		}
	}

	if(!is.null(hfile)){

		nFile=length(hfile)

		#Iteration among file
		cat("Output NUMERIC genotype...\n")
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
					if(i %% 1000 == 0)cat(paste("Number of Markers finished for theFile ", theFile, ": ", i, sep=""), "\n")
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
						GD = MRBLUP.Num(x=tt2[-c(1:11)], impute=SNP.impute)  
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
				myFile <- read.table(hfile[theFile], colClasses="character", sep=sep, head=FALSE, skip=1)
				nM <- nrow(myFile)
				write.table(myFile[, c(1, 3, 4)], paste(out, ".map", sep=""), append=TRUE, col.names=FALSE, row.names=FALSE, quote=FALSE, sep=sep)
				myFile <- myFile[, -c(1:11)];gc()
				myGDx <- apply(myFile, 1, function(x) MRBLUP.Num(x, impute=SNP.impute))
				myGDx <- t(myGDx)
				write.table(myGDx, paste(out, ".Numeric.txt", sep=""), append=TRUE, col.names=FALSE, row.names=FALSE, quote=FALSE, sep=sep)
				rm(myFile);rm(myGDx);gc()
				cat(paste("File: ", hfile[theFile], "; Total markers:", nM, " finished!\n", sep=""))
			}
		}
		cat("Preparation for NUMERIC data is done!\n")
	}

	#Transfer genotype data to .desc, .bin files
	if((!is.null(numfile))|(!is.null(hfile))){
		cat("Output BIG genotype...\n")
		if(is.null(numfile)){
			numfile <- paste(out, ".Numeric.txt", sep="")
			MAP <- read.table(paste(out, ".map", sep=""), head=FALSE, sep=sep)
		}else{
			MAP <- read.table(mapfile, head=FALSE, sep=sep)
			write.table(MAP, paste(out, ".map", sep=""), row.names=FALSE, col.names=FALSE, sep=sep, quote=FALSE)
		}
		nmarkers <- nrow(MAP); rm(MAP); gc()
		fileGenoCon <- file(description=numfile, open="r")
		tt2 <-readLines(fileGenoCon, n=1)
		tt2 <- unlist(strsplit(x=tt2, split=sep, fixed=TRUE))
		ns <- length(tt2)
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
				MRBLUP.Bar(i=i, n=nmarkers, fixed.points=FALSE)
				tt<-do.call(rbind, strsplit(x=tt, split=sep, fixed=TRUE))
				nn <- nrow(tt)
				#Identify end of file
				if(is.null(tt[1])) inGENOFile=FALSE
				if(inGENOFile){
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

			#Close GENO file
			close.connection(fileGenoCon)
		}
		if(priority == "speed"){
			myGeno <- read.big.matrix(numfile, type="char", head=FALSE, sep=sep,
				backingfile=bck, descriptorfile=dsc)
			rm("myGeno")
		}
		cat("Preparation for BIG data is done!\n")
	}
	gc()
	cat("MRBLUP data prepration accomplished successfully!\n")
}

MRBLUP.Num<-
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

MRBLUP.QTN.sel <-
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
	max.pos <-max(which(COR >= COR[1]))
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
			qtn.inc=qtn.new[[max.pos-1]]
			qtn.dec=NULL
		}else
		{
			qtn.dec=NULL
			
			#find the QTN that increase accuracy over the previous QTN and GBLUP
			#index <- ((cor.a-cor.new)<0) + (cor.a<cor.new[1])
			index <- ((cor.a-cor.new)<=0)
			dec.index <- which(index != 0)
			for(d in dec.index){
				if(d == 1){
					qtn.dec=c(qtn.dec, qtn.new[[d]])
				}else
				{
					qtn.dec=c(qtn.dec, setdiff(qtn.new[[d]], qtn.new[[d-1]]))
				}
			}
			qtn.inc=setdiff(qtn.new[[max.pos-1]], qtn.dec)
		}
	}
	return(list(qtn.inc=qtn.inc, qtn.dec=qtn.dec))
}

MRBLUP.CrossV <- 
function(
	phe, geno, K=NULL, CV=NULL, map=NULL, GWAS.model=NULL, vc.method="emma", Top.index=NULL, Top.perc=NULL, max.nQTN=TRUE, SNP.filter=0.5,
	cor.threshold=0.99, judge.threshold=0.9, sample.num=1, crv.num=5, cpu=1, theSeed=NULL, prior.QTN=NULL, K.method="center",
	bin.size=1000000, amplify=NULL, bisection=TRUE, bisection.loop=5, ref.gwas=FALSE, prior.model="QTN+K", GWAS.npc=3
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
# vc.method: method for variance components estimation("emma" or "gemmaU")								 #
# Top.index: a vector, a subset of top SNPs for each iteration are used as covariants					 #
# Top.perc: a vector, a subset of top SNPs for each iteration are amplified when calculating KINSHIP	 #
# max.nQTN: whether limits the max number of Top.index													 #
# SNP.filter: the P-values of T.test for SNPs which below this threshold will be deleted				 #
# cor.threshold: if the top SNP which is used as covariant is in high correlation with others,			 #
# 			it will be deleted																			 #
# sample.num: the sample number of cross-validation														 #
# crv.num: the cross number of cross-validation															 #
# judge.threshold: if the count of selected SNP for all iteration >= sample.num*crv.num*judge.threshold, #
# 			than it will be treated as covariance in final predicting model								 # 
# cpu: the number of CPU for calculation																 #
# theSeed: the random seed 																				 #
# prior.QTN: the prior QTNs which will be added as covariants, if provided prior QTNs, 	MRBLUP will not	 #
# 			optimize QTNs and model during cross-validation												 #
# prior.model: the prior MODEL(only "K", "QTN+K", "QTN" can be selected presently)						 #
# K.method: which algorithm will be applied("center", "scale", "vanraden")								 #
# bin.size: the size of each bin																		 #
# amplify: a vector, the base for LOG																	 #
# bisection: whether using bisection algorithm to optimize KINSHIP										 #
# bisection.loop: the max loop(iteration) number of bisection algorithm									 #
# ref.gwas: whether to do GWAS for reference population(if not, MRBLUP will merge all GWAS results of	 #
# 			cross-validation by mean)																	 #
# GWAS.npc=3: the number of PC that will be added as covariance to control population structure			 #
#--------------------------------------------------------------------------------------------------------#
	#make sure the type of system(windows, mac, linux)
	R.ver <- Sys.info()[['sysname']]
	wind <- R.ver == 'Windows'
	linux <- R.ver == 'Linux'
	mac <- (!linux) & (!wind)
    r.open <- !inherits(try(Revo.version,silent=TRUE),"try-error")
	if(is.null(cpu)){
		cat(" Please input the number of CPU:\n")
		cpu <- scan("", n=1, quiet=TRUE)
		cpu <- as.integer(cpu)
	}
	if(r.open &  mac & cpu > 1)	Sys.setenv("VECLIB_MAXIMUM_THREADS" = "1")
	if(wind)	cpu <- 1
	mkl <- try(getMKLthreads(), silent=TRUE)
	if(is.null(Top.index) & (is.null(Top.perc))) stop("One of 'Top.index' and 'Top.perc' must be setted")
	
	#make sure the type of R(base R, Open R), "checkpoint" package is installed in Open R by default
    #r.open <- "checkpoint" %in% rownames(installed.packages())
    r.open <- !inherits(try(Revo.version, silent=TRUE),"try-error")
	
	if(r.open &  mac & cpu > 1)	Sys.setenv("VECLIB_MAXIMUM_THREADS" = "1")
	
	N.Ind <- length(phe)
	NA.index <- which(is.na(phe))
	NA.Ind <- length(NA.index)
	Ind.index <- c(c(1:N.Ind)[-NA.index], NA.index)
	
	n.ref <- N.Ind-NA.Ind
	n.inf <- NA.Ind
	n.marker <- nrow(geno)
	binary <- length(unique(phe[-NA.index])) == 2

	if(!is.null(theSeed)) {set.seed(theSeed)}
	if(!is.null(prior.QTN))	Top.index <- NULL
	inf.index <- list()
	P.value <- list()
	if(!is.null(Top.index) | (!is.null(Top.perc) && (length(Top.perc) > 1 | length(amplify) > 1))){

		#get the inference population index for cross-validation
		for(s in 1:sample.num){
			CrossVindex <- sample(1:n.ref, n.ref)
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
			if(r.open)	try(setMKLthreads(MATH.cpu), silent=TRUE)
			ref.logic <- (i == sample.num * crv.num + 1)
			if(!ref.logic){
				gwas.ref.index <- Ind.index[c(1:n.ref)[-inf.index[[i]]]]
				P.ref <- phe
				P.ref[-gwas.ref.index] <- NA
			}else{
				P.ref <- phe
			}
			if(!ref.logic){
				pri <- paste(" GWAS of validations NO.", i, " finished ", sep="")
			}else{
				pri <- " GWAS of Total References finished "
			}
			if(GWAS.model == "GLM"){
				if(!is.null(SNP.filter)){
					GLM.gwas <- MRBLUP.GWAS(phe=P.ref, geno=geno[SNP.index, ], CV=CV, method="GLM", cpu=cpus, NPC=GWAS.npc, bar.head=pri, bar.tail="", bar.len=0)
					P.value <- vector("numeric", length=n.marker)
					P.value[SNP.index] <- GLM.gwas[, 2]
					P.value[-SNP.index] <- SNP.filter
				}else{
					GLM.gwas <- MRBLUP.GWAS(phe=P.ref, geno=geno, CV=CV,method="GLM", cpu=cpus, NPC=GWAS.npc, bar.head=pri, bar.tail="", bar.len=0)
					P.value <- GLM.gwas[, 2]
				}
				rm("GLM.gwas")
			}
			
			if(GWAS.model == "MLM"){
				if(!is.null(SNP.filter)){
					MLM.gwas <- MRBLUP.GWAS(phe=P.ref, geno=geno[SNP.index, ], CV=CV, K=K, method="MLM", vc.method="emma", cpu=cpus, NPC=GWAS.npc, bar.head=pri, bar.tail="", bar.len=0)
					P.value <- vector("numeric", length=n.marker)
					P.value[SNP.index] <- MLM.gwas[, 2]
					P.value[-SNP.index] <- SNP.filter
				}else{
					MLM.gwas <- MRBLUP.GWAS(phe=P.ref, geno=geno, CV=CV, K=K, method="MLM", vc.method="emma", cpu=cpus, NPC=GWAS.npc, bar.head=pri, bar.tail="", bar.len=0)
					P.value <- MLM.gwas[, 2]
					rm("MLM.gwas")
				}
			}
			rm(list=c("P.ref")); gc()
			P.value[is.na(P.value)] <- 1
			return(P.value)
		}
		
		GWAS.model.txt <- ifelse(GWAS.model=="GLM", "GLM(Generalized Linear Model)", "MLM(Mixed Linear Model)")
		cat(paste(" GWAS with model: ", GWAS.model.txt, "\n", sep=""))
		
		#the SNPs which are lower associated than 'SNP.filter' will be deleted to save time
		P.value.ref <- NULL
		if(!is.null(SNP.filter)){
			cat(paste(" Filtering SNPs with the threshold: ", SNP.filter, "\n", sep=""))
			GLM.gwas <- MRBLUP.GWAS(phe=phe, geno=geno, CV=CV, method="GLM", NPC=GWAS.npc, cpu=cpu, bar.head=" Filtering finished ", bar.tail="", bar.len=0)
			P.value.ref <- GLM.gwas[, 2]
			P.value.ref[is.na(P.value.ref)] <- 1
			SNP.index <- which(P.value.ref < SNP.filter)
			cat(paste(" Number of SNPs deleted: ", n.marker-length(SNP.index), "\n", sep=""))
			if(GWAS.model == "GLM" & ref.gwas){
				rm("GLM.gwas"); gc()
			}else{
				P.value.ref=NULL; rm("GLM.gwas"); gc()
			}
		}
		
		if(ref.gwas & is.null(P.value.ref)){
			gwas.num <- sample.num * crv.num + 1
		}else{
			gwas.num <- sample.num * crv.num
		}

		#because FarmCPU can not do parallel when testing SNP
		if(r.open)	MATH.cpu <- mkl

		if(wind & cpu > 1){
			cpus <- 1
			cat(" Multi-process of GWAS started...\n")
			max.cpu <- min(cpu, gwas.num)
			cl <- makeCluster(max.cpu, outfile = "Loop.log")
			registerDoParallel(cl)
			clusterExport(cl, varlist=c("MRBLUP.GWAS", "MRBLUP.EMMA.REML", "MRBLUP.Mix"))
			P.values <- foreach(x=1:gwas.num, .packages=c("bigmemory", "rfunctions")) %dopar% mult.run.gwas(x)
			stopCluster(cl)
			cat(" Multi-process done!\n")
		}else{
			cpus <- cpu
			P.values <- lapply(1:gwas.num, mult.run.gwas)
		}
		
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
		 # colnames(GWAS.res)[4] <- taxa
		 # wd <- getwd()
		 # setwd("../")
		 # CMplot(GWAS.res, plot.type="m",threshold=0.05,file="pdf",col=c("black","orange"))
		 # setwd(wd)
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
		cat(paste(" GWAS with model: ", GWAS.model, "...\n", sep=""))
		P.ref <- phe
		if(GWAS.model == "GLM"){
			GLM.gwas <- MRBLUP.GWAS(phe=P.ref, geno=geno, CV=CV, method="GLM", NPC=GWAS.npc, cpu=cpu, bar.head=" GWAS of References finished ", bar.tail="", bar.len=0)
			P.value.ref <- GLM.gwas[, 2]
			rm("GLM.gwas")
		}
		if(GWAS.model == "MLM"){
			MLM.gwas <- MRBLUP.GWAS(phe=P.ref, geno=geno, CV=CV, K=K, method="MLM", vc.method="emma", NPC=GWAS.npc, cpu=cpu, bar.head=" GWAS of References finished ", bar.tail="", bar.len=0)
			P.value.ref <- MLM.gwas[, 2]
			rm("MLM.gwas")
		}
		P.value.ref[is.na(P.value.ref)] <- 1
		rm(list=c("P.ref")); gc()
	}

	sam <- sample.num
	sample.num <- 1
	
	#Cross-validation to optimize pseudo QTN
	if(is.null(Top.index)){
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
		if(!0 %in% Top.index) {Top.index <- c(0, Top.index)}
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
		
		pseudo.QTN.k <-list()
		pseudo.QTN.qtn <-list()
		bin.qtn <- list()
		inc.QTN.store.k <- NULL
		inc.QTN.store.qtn <- NULL
		cor.pseudo.k <- NULL
		cor.pseudo.qtn <- NULL

		if(!is.null(bin.size)){
			cat(paste(" Set Genome into bins(", round(bin.size/1000000, 2), "Mb)...\n", sep=""))
			for(g in 1: length(P.value)){
				MRBLUP.gwas <- cbind(map, P.value[[g]])
				bin.qtn[[g]] <- MRBLUP.Bin(GWAS=MRBLUP.gwas, bin.size=bin.size, inclosure.size=max(Top.index), bin.num=1)
				
				#----------debug----------#
				 cat(" ", bin.qtn[[g]], "\n")
				#-------------------------#	
			}
			rm(MRBLUP.gwas);gc()
			cat(" Bins are confirmed!\n")
		}else{
			for(g in 1: length(P.value)){
				bin.qtn[[g]] <- order(P.value[[g]], decreasing=FALSE)[1:max(Top.index)]
			}
		}
		
		#if bin is too large, the number of bin.qtn would be less than the max of Top.index
		if(length(bin.qtn[[1]]) < max(Top.index)){
			cat(" (Warnings: bin size is too large, bin.qtns are not adequate for 'Top.index')\n")
			Top.index[Top.index > length(bin.qtn[[1]])] <- length(bin.qtn[[1]])
			Top.index <- sort(unique(Top.index))
		}
		cat(" Optimizing pseudo QTNs and MODEL...\n")
		#do parallel
		mult.run <- function(j, math.cpu=NULL){
		
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
			mult.ref.index <- Ind.index[c(1:n.ref)[-inf.index[[row.index]]]]
			mult.inf.index <- Ind.index[inf.index[[row.index]]]
			P.ref <- phe
			P.ref[-mult.ref.index] <- NA
			P.inf <- phe[mult.inf.index]

			if((Top.index[top]) == 0){
				gblup <- MRBLUP.Mix(phe=P.ref, K=K, CV=CV, vc.method=vc.method)
				gebv <- (CV %*% as.matrix(gblup$beta) + gblup$ebv)[mult.inf.index]
				if(binary){
					acc.k <- MRBLUP.stas.cal(P.inf, gebv, type="auc")
					acc.qtn <- MRBLUP.stas.cal(P.inf, gebv, type="auc")
				}else{
					acc.k <- MRBLUP.stas.cal(P.inf, gebv, type="cor")
					acc.qtn <- MRBLUP.stas.cal(P.inf, gebv, type="cor")
				}
				QTN.store <- NULL
				#cat(paste(" Cross-validation NO.", row.index, "; NQTN=", Top.index[top], "; Model='K'", "; ",mytext ,"=", round(acc.k, 4), paste(rep(" ", 5), collapse=""), "\r", sep=""))
			}else
			{
				myqtn <- bin.qtn[[row.index]][1:(Top.index[top])]
				if(length(myqtn) == 1){
					QTN.cv <- as.matrix(geno[myqtn, mult.ref.index])
				}else{
					QTN.cv <- t(geno[myqtn, mult.ref.index])
				}
				
				#remove the high correlated QTN(0.99)
				cor.cv <- cor(QTN.cv)
				cor.index <- abs(cor.cv) > cor.threshold
				b <- cor.cv * 0
				b[cor.index] <- 1
				bb <- 1-b
				bb[lower.tri(bb)] <- 1
				diag(bb) <- 1
				bd <- apply(bb, 2, prod)
				rm.index <- bd == 1
				myqtn <- myqtn[rm.index]
				QTN.cv <- QTN.cv[, rm.index]
				QTN.store <- myqtn
				if(length(myqtn) == 1){
					QTN.cv <- as.matrix(QTN.cv)
				}
				
				#calculate GLM model with QTNs
				qtn.eff <- MRBLUP.GLM(y=P.ref[mult.ref.index], X=CV[mult.ref.index, ], qtn.matrix=QTN.cv)$QTN.eff
				if(length(myqtn) == 1){
					gebv <- cbind(CV[mult.inf.index, ], geno[myqtn, mult.inf.index]) %*% as.matrix(qtn.eff)
					if(binary){
						acc.qtn <- MRBLUP.stas.cal(P.inf, gebv, type="auc")
					}else{
						acc.qtn <- MRBLUP.stas.cal(P.inf, gebv, type="cor")
					}
				}else{
					gebv <- cbind(CV[mult.inf.index, ], t(geno[myqtn, mult.inf.index])) %*% as.matrix(qtn.eff)
					if(binary){
						acc.qtn <- MRBLUP.stas.cal(P.inf, gebv, type="auc")
					}else{
						acc.qtn <- MRBLUP.stas.cal(P.inf, gebv, type="cor")
					}
				}
				
				#cat(paste(" Cross-validation NO.", row.index, "; NQTN=", Top.index[top], "; Model='QTN'", "; ", mytext, "=", round(acc.qtn, 4), paste(rep(" ", 5), collapse=""), "\r", sep=""))
				
				#calculate MLM model with QTNs
				MRBLUP.cv <- geno[myqtn, ]
				if(length(myqtn) != 1){
					MRBLUP.cv <- t(MRBLUP.cv)
				}
				MRBLUP.cv <- cbind(CV, as.matrix(MRBLUP.cv))
				gblup <- MRBLUP.Mix(phe=P.ref, vc.method=vc.method, CV=MRBLUP.cv, K=K)
				gebv <- (gblup$ebv + MRBLUP.cv %*% as.matrix(gblup$beta))[mult.inf.index]
				if(binary){
					acc.k <- MRBLUP.stas.cal(P.inf, gebv, type="auc")
				}else{
					acc.k <- MRBLUP.stas.cal(P.inf, gebv, type="cor")
				}
				#textx <- paste(" Cross-validation NO.", row.index, "; NQTN=", Top.index[top], "; Model='QTN+K'", "; ", mytext, "=", round(acc.k, 4), sep="")
				#cat( textx, paste(rep(" ", 5), collapse=""), "\r")
			}
			print.f(j)
			rm(list=c("P.ref", "P.inf", "gblup", "gebv")); gc()
			return(list(acc.k=acc.k, acc.qtn=acc.qtn, QTN.store=QTN.store))
		}
		iterationN <- ((sample.num * crv.num) * length(Top.index))
		cat(paste(" Total iteration number:", iterationN))
		if(cpu == 1){
			print.f <- function(i){MRBLUP.Bar(i=i, n=iterationN, type="type1", symbol.len=0, symbol.head="Cross-validation Finished_", symbol.tail="")}
			mult.res <- lapply(1:iterationN, mult.run); gc()
		}else
		{
			#1：mclapply doesn't work at windows system
			#2: there are always some errors when use multi-process in microsoft r open
			
			#get the user-specific cpu number
 			if(wind){
				print.f <- function(i){MRBLUP.Bar(i=i, n=iterationN, type="type1", symbol.len=0, symbol.head="Cross-validation Finished_", symbol.tail="")}
				#print(paste("Cross-validation NO.", paste(c(1:(sample.num * crv.num)), collapse=", "), sep=""))
				cat(" Multi-process started...\n")
				cat(" (Note: There needs to wait some time! See iteration details in 'Loop.log')\n")
				cl <- makeCluster(cpu, outfile = "Loop.log")
				registerDoParallel(cl)
				clusterExport(cl, varlist=c("MRBLUP.Mix", "MRBLUP.EMMA.REML", "MRBLUP.GEMMA.REML", "MRBLUP.GEMMA.VC",
					"MRBLUP.stas.cal", "MRBLUP.Kin", "MRBLUP.QTN.rm", "MRBLUP.GLM"))
				mult.res <- foreach(x=1:iterationN,
                .packages=c("bigmemory", "rfunctions")) %dopar% mult.run(x, math.cpu=MATH.cpu)
				if(is.null(Top.perc))	stopCluster(cl)
				cat(" Multi-process done!\n")
			}else{
				tmpf.name <- tempfile()
				tmpf <- fifo(tmpf.name, open="w+b", blocking=TRUE)
				
				# print.f <- function(i){writeBin(1, tmpf)}
				# MRBLUP.Bar(n=iterationN, type="type2", tmp.file=tmpf, fixed.points=FALSE, symbol.len=0, symbol.head=" Cross-validation Finished_", symbol.tail="")
				
				writeBin(0, tmpf)
				print.f <- function(i){MRBLUP.Bar(n=iterationN, type="type3", tmp.file=tmpf, fixed.points=FALSE, symbol.len=0, symbol.head=" Cross-validation Finished_", symbol.tail="")}
				
				if(r.open)	setMKLthreads(1)
				mult.res <- parallel::mclapply(1:iterationN, mult.run, mc.cores = cpu)
				if(r.open)	setMKLthreads(mkl)
				close(tmpf); unlink(tmpf.name)
			}
			gc()
		}
		res <- do.call(rbind, mult.res)
		cor.store.k <- matrix(unlist(res[, 1]), length(Top.index))
		
		cor.store.qtn <- matrix(unlist(res[, 2]), length(Top.index))
		qtn.listT <- res[, 3]
		
		#Select the qtns which can increase the accuracy for each validation
		for(scv in 1:(sample.num * crv.num)){
			mylist <- qtn.listT[(length(Top.index) * (scv-1) + 2):(length(Top.index) * scv)]
			mrblup.qtn <- MRBLUP.QTN.sel(COR=cor.store.k[, scv], QTN.list=mylist)$qtn.inc
			
			#----------debug----------#
			#print(mrblup.qtn)
			#-------------------------#	
			
			inc.QTN.store.k<- c(inc.QTN.store.k, mrblup.qtn)
			inc.QTN.store.qtn<- c(inc.QTN.store.qtn, MRBLUP.QTN.sel(COR=cor.store.qtn[, scv], QTN.list=mylist)$qtn.inc)
		}
		
		#merge the selected QTNs by bins(only one QTN is permitted in each bin)
		if(!is.null(inc.QTN.store.k)){
		
			#----------debug----------#
			#print(table(inc.QTN.store.k))
			#-------------------------#
			
			count.inc.QTN.store.k <- table(inc.QTN.store.k)
			inc.QTN.store.k <- as.numeric(names(count.inc.QTN.store.k))
			
			#choose the replacement 'method' for each effective QTN
			
			#get the P-value from merge the GWAS of cross-validation:
			#inc.QTN.store.k.p <- do.call(rbind, lapply(P.value, function(x) x[c(inc.QTN.store.k)]))
			#inc.QTN.store.k.p <- apply(inc.QTN.store.k.p, 2, P.select)
			
			#get the P-value from t-test
			# inc.QTN.store.k.p <- unlist(lapply(1: length(inc.QTN.store.k), function(x){
				# P <- cor.test(refphe[,2],refgeno[inc.QTN.store.k[x],])$p.value
				# return(P)
			# }))
			
			#get the P-value from selected model
			inc.QTN.store.k.p <- P.value.ref[inc.QTN.store.k]
			
			#order the QTNs by p-value
			ordered.inc.QTN.store.k <- order(inc.QTN.store.k.p, decreasing=FALSE)
			
			inc.QTN.store.k <- inc.QTN.store.k[ordered.inc.QTN.store.k]
			count.inc.QTN.store.k <- count.inc.QTN.store.k[ordered.inc.QTN.store.k]
			inc.QTN.store.k.map <- map[inc.QTN.store.k,]
			inc.QTN.store.k.x <- NULL
			count.inc.QTN.store.k.x <- NULL
			for(INC in 1: length(inc.QTN.store.k)){
				
				#if inc.qtn.store.k is null or length one, stop loop!
				if(length(inc.QTN.store.k) %in% c(0,1)) {
					inc.QTN.store.k <- c(inc.QTN.store.k.x, inc.QTN.store.k)
					count.inc.QTN.store.k <- c(count.inc.QTN.store.k.x, count.inc.QTN.store.k)
					break()
				}
				inc.QTN.store.k.x <- c(inc.QTN.store.k.x, inc.QTN.store.k[1])
				index1 <- inc.QTN.store.k.map[,2] == inc.QTN.store.k.map[1,2]
				index2 <- (inc.QTN.store.k.map[, 3] >= (inc.QTN.store.k.map[1, 3]-bin.size/2)) & (inc.QTN.store.k.map[, 3] <= (inc.QTN.store.k.map[1, 3] + bin.size/2))
				count.inc.QTN.store.k.x <- c(count.inc.QTN.store.k.x, sum(count.inc.QTN.store.k[which(index1 & index2)]))
				
				#delete the selected qtns
				inc.QTN.store.k <- inc.QTN.store.k[!(index1 & index2)]
				inc.QTN.store.k.map <- inc.QTN.store.k.map[!(index1 & index2), ]
				count.inc.QTN.store.k <- count.inc.QTN.store.k[!(index1 & index2)]
				names(count.inc.QTN.store.k) <- NULL
			}
			
			#----------------------debug----------------------#
			#print(rbind(inc.QTN.store.k, count.inc.QTN.store.k))
			#-------------------------------------------------#
			
			inc.QTN.store.k <- inc.QTN.store.k[which(count.inc.QTN.store.k >= (floor(sample.num * crv.num * judge.threshold)))]
			
			if(length(inc.QTN.store.k) == 0)	inc.QTN.store.k=NULL
		}
		
		#-----------debug-----------#
		#print(inc.QTN.store.k)
		#---------------------------#
		
		if(!is.null(inc.QTN.store.qtn)){
			count.inc.QTN.store.qtn <- table(inc.QTN.store.qtn)
			inc.QTN.store.qtn <- as.numeric(names(count.inc.QTN.store.qtn))
			# inc.QTN.store.qtn.p <- do.call(rbind, lapply(P.value, function(x) x[c(inc.QTN.store.qtn)]))
			# inc.QTN.store.qtn.p <- apply(inc.QTN.store.qtn.p, 2, P.select)
			# inc.QTN.store.qtn.p <- unlist(lapply(1: length(inc.QTN.store.qtn), function(x){
				# P <- cor.test(refphe[,2],refgeno[inc.QTN.store.qtn[x],])$p.value
				# return(P)
			# }))
			inc.QTN.store.qtn.p <- P.value.ref[inc.QTN.store.qtn]
			ordered.inc.QTN.store.qtn <- order(inc.QTN.store.qtn.p, decreasing=FALSE)
			inc.QTN.store.qtn <- inc.QTN.store.qtn[ordered.inc.QTN.store.qtn]
			count.inc.QTN.store.qtn <- count.inc.QTN.store.qtn[ordered.inc.QTN.store.qtn]
			inc.QTN.store.qtn.map <- map[inc.QTN.store.qtn,]
			inc.QTN.store.qtn.x <- NULL
			count.inc.QTN.store.qtn.x <- NULL
			for(INC in 1: length(inc.QTN.store.qtn)){
				if(length(inc.QTN.store.qtn) %in% c(0,1)) {
					inc.QTN.store.qtn <- c(inc.QTN.store.qtn.x, inc.QTN.store.qtn)
					count.inc.QTN.store.qtn <- c(count.inc.QTN.store.qtn.x, count.inc.QTN.store.qtn)
					break()
				}
				inc.QTN.store.qtn.x <- c(inc.QTN.store.qtn.x, inc.QTN.store.qtn[1])
				index1 <- inc.QTN.store.qtn.map[,2] == inc.QTN.store.qtn.map[1,2]
				index2 <- (inc.QTN.store.qtn.map[, 3] >= (inc.QTN.store.qtn.map[1, 3]-bin.size/2)) & (inc.QTN.store.qtn.map[, 3] <= (inc.QTN.store.qtn.map[1, 3] + bin.size/2))
				count.inc.QTN.store.qtn.x <- c(count.inc.QTN.store.qtn.x, sum(count.inc.QTN.store.qtn[which(index1 & index2)]))
				inc.QTN.store.qtn <- inc.QTN.store.qtn[!(index1 & index2)]
				inc.QTN.store.qtn.map <- inc.QTN.store.qtn.map[!(index1 & index2), ]
				count.inc.QTN.store.qtn <- count.inc.QTN.store.qtn[!(index1 & index2)]
			}
			inc.QTN.store.qtn <- inc.QTN.store.qtn[which(count.inc.QTN.store.qtn >= (floor(sample.num * crv.num * judge.threshold)))]
			if(length(inc.QTN.store.qtn) == 0)	inc.QTN.store.qtn=NULL
		}
		
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
		cat("\r", "The optimized MODEL:", paste(rep(" ", 20), collapse=""), "\n")
		cat(" ", cross.model, "\n")
		cat(" The optimized QTNs:\n")
		if(is.null(cross.QTN)){
			cat(" NULL", "\n")
		}else{
			cat(" ", cross.QTN, "\n")
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
			mult.ref.index <- Ind.index[c(1:n.ref)[-inf.index[[cv]]]]
			mult.inf.index <- Ind.index[inf.index[[cv]]]
			P.ref <- phe
			P.ref[-mult.ref.index] <- NA
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
				max.wt <- max(top.wt)
				
				#set the weight of cross.QTN to 0
				#top.wt[cross.QTN] <- 0
				
				#-----------debug-----------#
                #print(range(top.wt))
                #print(sum(is.na(top.wt)))
				#---------------------------#

				Kt <- MRBLUP.Kin(geno[mytop.perc, ], weight=top.wt, SUM=SUM, type=K.method)
				Kt <- K + Kt
				
				#-----------debug-----------#
				# print(Kt[1:3,1:3])
				# if(j <= 20)	write.csv(myK, paste(j,"_",amplify[amp.index], "_err.K.csv", sep=""), row.names=FALSE)
				#---------------------------#

				rm(list=c("mytop.order", "mytop.perc", "top.wt")); gc()
			}
			if(!is.null(Kt)){
				if(!is.null(cross.QTN)){
					if(length(cross.QTN) == 1){
						QTN.cv <- as.matrix(geno[cross.QTN, ])
					}else{
						QTN.cv <- t(geno[cross.QTN, ])
					}
				}else
				{
					QTN.cv <- NULL
				}
				MRBLUP.cv <- cbind(CV, QTN.cv)
				gblup <- MRBLUP.Mix(phe=P.ref, vc.method=vc.method, CV=MRBLUP.cv, K=Kt)
				MRBLUP.ll <- gblup$LL
				gebv <- (gblup$ebv + MRBLUP.cv %*% as.matrix(gblup$beta))[mult.inf.index]
				if(binary){
					acc<- MRBLUP.stas.cal(P.inf, gebv, type="auc")
				}else{
					acc<- MRBLUP.stas.cal(P.inf, gebv, type="cor")
				}
				rm(list=c("P.ref", "P.inf", "Kt", "gblup", "gebv")); gc()
				#cat(paste(" Cross-validation NO.", cv, "; QTN=", ceiling(n.marker * Top.perc[top]), "(", Top.perc[top] * 100,
                #"%); Logx=", round(amplify[amp.index], 4),"; ",mytext , "=", round(acc, 4), paste(rep(" ", 5), collapse=""), "\r", sep=""))
			}else{
				acc <- NA; max.wt <- NA; MRBLUP.ll <- NA
			}
			print.f(j)
			return(list(acc=acc,max.wt=max.wt,MRBLUP.ll=MRBLUP.ll))
		}
		iterationN <- ((sample.num * crv.num) * length(Top.perc) * length(amplify))
		cat(paste(" Total iteration number:", iterationN))
		MATH.cpu <- mkl
		if(cpu == 1){
			print.f <- function(i){MRBLUP.Bar(i=i, n=iterationN, type="type1", symbol.len=0, symbol.head="Cross-validation Finished_", symbol.tail="")}
			mult.res <- lapply(1 : iterationN, mult.run)
		}else
		{
			if(wind){
				print.f <- function(i){MRBLUP.Bar(i=i, n=iterationN, type="type1", symbol.len=0, symbol.head="Cross-validation Finished_", symbol.tail="")}
				cat(" Multi-process started...\n")
				cat(" (Note: There needs to wait some time! See iteration details in 'Loop.log')\n")
				if(is.null(Top.index)){
					cl <- makeCluster(cpu, outfile = "Loop.log")
					registerDoParallel(cl)
					clusterExport(cl, varlist=c("MRBLUP.Mix", "MRBLUP.EMMA.REML", "MRBLUP.stas.cal", "MRBLUP.Kin", "MRBLUP.QTN.rm", "MRBLUP.GLM"))
					#Exp.packages <- clusterEvalQ(cl, c(library(bigmemory), library(rfunctions)))
				}
				mult.res <- foreach(x=1:iterationN,
                .packages=c("bigmemory", "rfunctions")) %dopar% mult.run(x, math.cpu=MATH.cpu)
				if(!bisection)	stopCluster(cl)
				cat(" Multi-process done!\n")
			}else
			{
				tmpf.name <- tempfile()
				tmpf <- fifo(tmpf.name, open="w+b", blocking=TRUE)
				
				# print.f <- function(i){writeBin(1, tmpf)}
				# MRBLUP.Bar(n=iterationN, type="type2", tmp.file=tmpf, fixed.points=FALSE, symbol.len=0, symbol.head=" Cross-validation Finished_", symbol.tail="")
				
				writeBin(0, tmpf)
				print.f <- function(i){MRBLUP.Bar(n=iterationN, type="type3", tmp.file=tmpf, fixed.points=FALSE, symbol.len=0, symbol.head=" Cross-validation Finished_", symbol.tail="")}
				
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
		if(bisection){
			cat(" Bisection algorithm started...\n")
			TOP <- Top.perc
			AMP <- amplify
			K.COR <- K.cor.store
			for(loop in 1 : bisection.loop){
				K.cor.mean <- apply(K.cor.store, 1, mean)
				#K.cor.mean <- apply(K.cor.store, 1, median)
				top.cor <- K.cor.store[which.max(K.cor.mean), ]
				
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
				MATH.cpu <- mkl
				cat("\r", paste("Loop ", loop, " of ", bisection.loop, " selected interval: [", round(top1 * 100, 4), "%,", round(top2 * 100, 4), "%]; ", "[", round(amp1, 4), ",", round(amp2, 4), "]\n", sep=""))
				cat(paste(" Total iteration number:", iterationN))
				if(cpu == 1){
					print.f <- function(i){MRBLUP.Bar(i=i, n=iterationN, type="type1", symbol.len=0, symbol.head="Cross-validation Finished_", symbol.tail="")}
					mult.res <- lapply(1:iterationN, mult.run)
				}else
				{
					if(wind){
						print.f <- function(i){MRBLUP.Bar(i=i, n=iterationN, type="type1", symbol.len=0, symbol.head="Cross-validation Finished_", symbol.tail="")}
						mult.res <- foreach(x=1:iterationN,
						.packages=c("bigmemory", "rfunctions")) %dopar% mult.run(x, math.cpu=MATH.cpu)
						if(loop == bisection.loop)	stopCluster(cl)
					}else
					{
						tmpf.name <- tempfile()
						tmpf <- fifo(tmpf.name, open="w+b", blocking=TRUE)
						
						# print.f <- function(i){writeBin(1, tmpf)}
						# MRBLUP.Bar(n=iterationN, type="type2", tmp.file=tmpf, fixed.points=FALSE, symbol.len=0, symbol.head=" Cross-validation Finished_", symbol.tail="")
						
						writeBin(0, tmpf)
						print.f <- function(i){MRBLUP.Bar(n=iterationN, type="type3", tmp.file=tmpf, fixed.points=FALSE, symbol.len=0, symbol.head=" Cross-validation Finished_", symbol.tail="")}
				
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
				if(length(amplify) == 1){top.index <- c(top, Top.perc)}else if(length(Top.perc) == 1)
					{top.index <- c(top, rep(Top.perc, 2))}else{top.index <- c(top, rep(Top.perc, c(2, 2)))}
				top.amplify <- c(amp, rep(amplify, length(Top.perc)))
				colnames(K.cor.store) <- 1:(sample.num * crv.num)
				rownames(K.cor.store) <- paste(top.index, top.amplify,sep="_")
				if(length(Top.perc) != 1)	Top.perc <- c(top, Top.perc)
				if(length(amplify) != 1)	amplify <- c(amp, amplify)
				
				#-----------debug-----------#
				#print(K.cor.store)
				#---------------------------#

				rm(list=c("mult.store", "mult.res")); gc()
			}
			K.cor.mean <- apply(K.cor.store, 1, mean)
			K.cor.median <- apply(K.cor.store, 1, median)
			top.index.max <- top.index[which.max(K.cor.mean)]
			K.cor.median.max <- K.cor.median[which.max(K.cor.mean)]
			amp.max <- top.amplify[which.max(K.cor.mean)]
		}
		
        if(
			(max(K.cor.mean) > mean(cor_qtn_k)) & (K.cor.median.max >= median(cor_qtn_k))
			#& 
			#sum((K.cor.store[which.max(K.cor.mean), ] - cor_qtn_k) > 0) >= floor((sample.num * crv.num) * judge.threshold)
		)
		{
			K.opt <- TRUE
			cat("\r", paste("Posterior top percentage: ", paste(rep(" ", 20), collapse=""), "\n", sep=""))
			cat(" ", paste(top.index.max * 100, "%", sep=""), "\n")
			cat(paste(" Posterior Logx: ", "\n", sep=""))
			cat(" ", round(amp.max, 4), "\n")
		}else{
			K.opt <- FALSE
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
			
			#top.wt <- top.wt*1.5
			
			#-----------debug-----------#
            # print(range(top.wt))
			#---------------------------#
			
			cat(" Calculating optimized KINSHIP...\n")
			optK <- MRBLUP.Kin(geno[mytop.perc,], type=K.method, weight = top.wt, SUM=SUM)
            K <- K + optK
            
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
			optK <- MRBLUP.Kin(geno[mytop.perc,], type="center", weight = top.wt, SUM=nrow(geno))
            K <- K + optK
			rm(list=c("mytop.order", "mytop.perc", "top.wt", "optK"))
			rm(list=c("P.value.ref")); gc()
		}else{
			cat(" No Kinship optimization\n")
		}
	}
	cat(" Cross-validation DONE!\n")
	return(list(cross.QTN = cross.QTN, cross.model=cross.model, cross.k=K))
}

MRBLUP.GLM <- 
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
	XX <- try(crossprodcpp(X), silent=TRUE)
	if(inherits(XX, "try-error")){
		warning("Package 'rfunctions' is not installed!")
		XX <- crossprod(X)
	}
	XY<-crossprod(X, Y)
	YY <- try(crossprodcpp(Y), silent=TRUE)
	if(inherits(YY, "try-error")){
		warning("Package 'rfunctions' is not installed!")
		YY <- crossprod(Y)
	}
	XXi <- try(solve(XX), silent = TRUE)
	if(inherits(XXi, "try-error")){
		if(r.open){
			XXi <- MASS::ginv(XX)
		}else{
			XXi <- try(geninv(XX), silent = TRUE)
			if(inherits(XXi, "try-error")){
				warning("Package 'rfunctions' is not installed!")
				XXi <- MASS::ginv(XX)
			}
		}
	}
	beta<-crossprod(XXi, XY)
	QTN.eff<-as.numeric(beta)
	return(list(QTN.eff=QTN.eff))
}

MRBLUP.GWAS <- 
function(
	phe, geno, K=NULL, CV=NULL, NPC=NULL, REML=NULL, cpu=NULL, vc.method="emma", method="MLM", bar.head="|", bar.tail=">", bar.len=50
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
# vc.method: method for variance components estimation("emma" or "gemmaU")								 #
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
	if(r.open)	math.cpu <- try(getMKLthreads(), silent=TRUE)
	if(is.null(cpu)){
		cat(" Please input the number of CPU:\n")
		cpu <- scan("", n=1, quiet=TRUE)
		cpu <- as.integer(cpu)
	}
	if(r.open &  mac & cpu > 1)	Sys.setenv("VECLIB_MAXIMUM_THREADS" = "1")
	if(wind)	cpu <- 1

	Ind.index <- which(!is.na(phe))
	ys <- as.numeric(phe[Ind.index])
	if(length(Ind.index) == length(phe)){
		geno <- as.matrix(geno)
	}else{
		geno <- geno[, Ind.index]
	}
	n <- ncol(geno)
	m <- nrow(geno)
	if(method == "MLM"){
		if(is.null(K)){
			K <- MRBLUP.Kin(geno)
		}else{
			K <- K[Ind.index, Ind.index]
		}
	}

	if(!is.null(NPC)){
		pca <- prcomp((t(geno)))$x[, 1:NPC]
	}else{
		pca <- NULL
	}
	if(is.null(CV)){
		X0 <- cbind(matrix(1, n), pca)
	}else{
		X0 <- cbind(pca, CV[Ind.index, ])
	}
	if(is.null(REML) & method == "MLM"){
		if(vc.method == "gemmaU") REML <- MRBLUP.GEMMA.REML(ys, X=X0, K=K, root=FALSE)
		if(vc.method == "emma") REML <- MRBLUP.EMMA.REML(ys, X=X0, K=K)
	}
	q0 <- ncol(X0)
	iXX <- matrix(NA,q0+1,q0+1)
	Xt <- matrix(NA,n, q0+1)
	
	#parallel function for MLM model
	eff.mlm <- function(i){
		SNP <- geno[i, ]
		#if(min(SNP) == max(SNP)){
			#effect <- 0;p <- 1
		#}else
		#{
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
			print.f(i)
		#}
		return(list(effect = effect, p = p))
	}
	
	#parallel function for GLM model
	eff.glm <- function(i){
		SNP <- geno[i, ]
		#if(min(SNP) == max(SNP)){
			#effect <- 0;p <- 1
		#}else
		#{
			#Process the edge (marker effects)
			sy <- crossprod(SNP,y)
			ss <- crossprod(SNP)
			xs <- crossprod(X0,SNP)
			
			B21 <- crossprod(xs, X0X0i)
			t2 <- B21 %*% xs
			B22 <- ss - t2
			invB22 <- 1/B22
			NeginvB22B21 <- crossprod(-invB22,B21)
			B21B21 <- crossprod(B21)
			iXX11 <- X0X0i + as.numeric(invB22) * B21B21
			
			#Derive inverse of LHS with partationed matrix
			iXX[1:q0,1:q0] <- iXX11
			iXX[(q0+1),(q0+1)] <- invB22
			iXX[(q0+1),1:q0] <- NeginvB22B21
			iXX[1:q0,(q0+1)] <- NeginvB22B21
			df <- n-q0-1
			rhs <- c(X0Y,sy)
			effect <- crossprod(iXX,rhs)
			ve <- (YY-crossprod(effect,rhs))/df #this is a scaler
			effect <- effect[q0+1]
			t.value <- effect/sqrt(iXX[q0+1, q0+1] * ve)
			p <- 2 * pt(abs(t.value), df, lower.tail=FALSE)
			print.f(i)
		#}
		return(list(effect=effect, p=p))
    }
    
    if(method == "MLM"){
        ves <- REML$ve
        vgs <- REML$vg
        lambda <- ves/vgs
        if(ncol(K) < 2000) {
            try(setMKLthreads(1), silent=TRUE)
        }
        eig <- eigen(K, symmetric=TRUE)
        if(ncol(K) < 2000) {
            try(setMKLthreads(math.cpu), silent=TRUE)
        }
        U <- eig$vectors * matrix(sqrt(1/(eig$values + lambda)), n, length(eig$values), byrow=TRUE)
		y <- matrix(ys)
		yt <- crossprod(U, y)
		X0t <- crossprod(U, X0)
		X0X0 <- crossprod(X0t)
		X0Y <- crossprod(X0t,yt)
		iX0X0 <- solve(X0X0)
		Xt[1:n,1:q0] <- X0t
    }
	
    if(method == "GLM"){
		y <- matrix(ys)
        X0X0 <- crossprod(X0)
		X0Y <- crossprod(X0,y)
		YY <- crossprod(y)
		X0X0i <- solve(X0X0)
    }
    if(cpu == 1){
		print.f <- function(i){MRBLUP.Bar(i=i, n=m, type="type1", symbol.len=bar.len, symbol.head=bar.head, symbol.tail=bar.tail)}
		if(method == "MLM") results <- lapply(1:m, eff.mlm)
		if(method == "GLM") results <- lapply(1:m, eff.glm)
	}else
	{
		if(method == "MLM"){
			if(wind){
				print.f <- function(i){MRBLUP.Bar(i=i, n=m, type="type1", symbol.len=bar.len, symbol.head=bar.head, symbol.tail=bar.tail)}
				#foreach function
				#print(" *  *  *  *  *  * test foreach parallel *  *  *  *  *  * ")
				#cl <- makeCluster(cpu)
				#registerDoParallel(cl)
				#Exp.packages <- clusterEvalQ(cl, c(library(bigmemory)))
				#results <- foreach(x=1:m) %dopar% eff.mlm(x)
				#stopCluster(cl)
				
				#print(" *  *  *  *  *  * test parLapply parallel *  *  *  *  *  * ")
				cl <- makeCluster(getOption("cl.cores", cpu))
				clusterExport(cl, varlist=c("geno", "yt", "X0", "U", "vgs", "ves", "math.cpu"), envir=environment())
				Exp.packages <- clusterEvalQ(cl, c(library(bigmemory),library(rfunctions)))
				results <- parallel::parLapply(cl, 1:m, eff.mlm)
				stopCluster(cl)
			}else{
				tmpf.name <- tempfile()
				tmpf <- fifo(tmpf.name, open="w+b", blocking=TRUE)
				
				# print.f <- function(i){writeBin(1, tmpf)}
				# MRBLUP.Bar(n=m, type="type2", tmp.file=tmpf, symbol.len=bar.len, symbol.head=bar.head, symbol.tail=bar.tail)
				
				writeBin(0, tmpf)
				print.f <- function(i){MRBLUP.Bar(n=m, type="type3", tmp.file=tmpf, symbol.len=bar.len, symbol.head=bar.head, symbol.tail=bar.tail)}
				
				if(r.open){
					if(R.ver == 'Linux' && r.open)	try(setMKLthreads(1), silent=TRUE)
					results <- parallel::mclapply(1:m, eff.mlm, mc.cores=cpu)
					if(R.ver == 'Linux' && r.open)	try(setMKLthreads(math.cpu), silent=TRUE)
				}else{
					results <- parallel::mclapply(1:m, eff.mlm, mc.cores=cpu)
				}
				#Sys.sleep(0.5)
				close(tmpf); unlink(tmpf.name); cat('\n');
			}
		}
		if(method  == "GLM"){
			if(wind){
				print.f <- function(i){MRBLUP.Bar(i=i, n=m, type="type1", symbol.len=bar.len, symbol.head=bar.head, symbol.tail=bar.tail)}
				#parLapply function
				#print(" *  *  *  *  *  * test parLapply parallel *  *  *  *  *  * ")
				cl <- makeCluster(getOption("cl.cores", cpu))
				clusterExport(cl, varlist=c("geno", "ys", "X0", "math.cpu"), envir=environment())
				Exp.packages <- clusterEvalQ(cl, c(library(bigmemory),library(rfunctions)))
				results <- parallel::parLapply(cl, 1:m, eff.glm)
				stopCluster(cl)
			}else{
				tmpf.name <- tempfile()
				tmpf <- fifo(tmpf.name, open="w+b", blocking=T)
				
				# print.f <- function(i){writeBin(1, tmpf)}
				# MRBLUP.Bar(n=m, type="type2", tmp.file=tmpf, symbol.len=bar.len, symbol.head=bar.head, symbol.tail=bar.tail)

				writeBin(0, tmpf)
				print.f <- function(i){MRBLUP.Bar(n=m, type="type3", tmp.file=tmpf, symbol.len=bar.len, symbol.head=bar.head, symbol.tail=bar.tail)}
				
				if(r.open){
					if(R.ver == 'Linux' && r.open)	try(setMKLthreads(1), silent=TRUE)
					results <- parallel::mclapply(1:m, eff.glm, mc.cores=cpu)
					if(R.ver == 'Linux' && r.open)	try(setMKLthreads(math.cpu), silent=TRUE)
				}else{
					results <- parallel::mclapply(1:m, eff.glm, mc.cores=cpu)
				}
				#Sys.sleep(0.5)
				close(tmpf); unlink(tmpf.name); cat('\n');
			}
		}
	}
	if(is.list(results)) results <- matrix(unlist(results), m, byrow=TRUE)
	return(results)
}

MRBLUP.copy <- 
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
			memo <- "MRBLUP.temp"
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

MRBLUP <- 
function(
	pfile="", gfile="", kfile=NULL, cfile=NULL, pheno=1, GWAS.model="MLM", GWAS.npc=NULL, prior.QTN=NULL, prior.model="QTN+K", vc.method="emma", 
	K.method="center",Top.perc=c(1e-4, 1e-3, 1e-2, 1e-1), Top.index=1:20, Logx=c(1.01, 1.11, exp(1), 10),
	bin.size=1000000, max.nQTN=TRUE, sample.num=4, SNP.filter=NULL, crv.num=5, cor.threshold=0.3, judge.threshold=0.9, bisection=TRUE,
	bisection.loop=8, ref.gwas=FALSE, theSeed=666666, file.output=TRUE, cpu=NULL
)
{
#--------------------------------------------------------------------------------------------------------#
# Object: To perform GP/GS with MRBLUP(Genomic Prediction/Selection)									 #
#																										 #
# Input:																								 #
# pfile: phenotype file, one column for a trait, the name of each column must be provided(NA is allowed) #
# gfile: genotype files, including "gfile.geno.desc", "gfile.geno.bin" and "gfile.map"					 #
# kfile: n*n, optional, provided KINSHIP file for all individuals										 #
# cfile: n*x, optional, the provided covariance file													 #
# pheno：specify phenotype column in the phenotype file(default 1)										 #
# GWAS.model: which model will be used for GWAS(only "GLM" and "MLM" can be selected presently)			 #
# GWAS.npc: the number of PC that will be added as covariance to control population structure			 #
# prior.QTN: the prior QTNs which will be added as covariants, if provided prior QTNs, 	MRBLUP will not	 #
# 			optimize QTNs and model during cross-validation												 #
# prior.model: the prior Model for the prior.QTN that added as covariants								 #
# vc.method: method for variance components estimation("emma" or "gemmaU")								 #
# K.method: which algorithm will be applied("center", "scale", "vanraden")								 #
# Top.perc: a vector, a subset of top SNPs for each iteration are amplified when calculating KINSHIP	 #
# Top.index: a vector, a subset of top SNPs for each iteration are used as covariants					 #
# Logx: a vector, the base for LOG																		 #
# bin.size: the size of each bin																		 #
# max.nQTN: whether limits the max number of Top.index													 #
# sample.num: the sample number of cross-validation														 #
# crv.num: the cross number of cross-validation 														 #
# SNP.filter: the SNPs whose P-value below this threshold will be deleted								 #
# cor.threshold: if the top SNP which is used as covariant is in high correlation with others,			 #
# 			it will be deleted																			 #
# judge.threshold: if the count of selected SNP for all iteration >= sample.num*crv.num*judge.threshold, #
# 			than it will be treated as covariance in final predicting model								 # 
# bisection: whether using bisection algorithm to optimize KINSHIP										 #
# bisection.loop: the max loop(iteration) number of bisection algorithm									 #
# ref.gwas: whether to do GWAS for reference population(if not, MRBLUP will merge all GWAS results of	 #
# 			cross-validation by mean)																	 #
# theSeed: the random seed																				 #
# file.output: whether to write the predicted values in file											 #
# cpu: the number of CPU for calculation																 #
#--------------------------------------------------------------------------------------------------------#
	#print the version of MRBLUP
	MRBLUP.version()
    time1 <- as.numeric(Sys.time())

    #check the parameters
    if(!GWAS.model %in% c("GLM", "MLM")) stop("Please select the right GWAS model!('GLM', 'MLM')")
	if(!K.method %in% c("scale", "vanraden", "center")) stop("Please select the right K.method!('scale', 'vanraden', 'center')")
	if(!vc.method %in% c("emma", "gemmaU"))   stop("Please choose only one of c('emma', 'gemmaU')!")
    if(crv.num < 2) stop("'crv.num' must be bigger than 2!")
    if(sample.num < 1) stop("'sample.num' must be bigger than 1!")
	if(is.null(cpu)){
		cat(" Please input the number of CPU:\n")
		cpu <- scan("", n=1, quiet=TRUE)
		cpu <- as.integer(cpu)
	}
	
	R.ver <- Sys.info()[['sysname']]
	wind <- R.ver == 'Windows'
	linux <- R.ver == 'Linux'
	mac <- (!linux) & (!wind)
    r.open <- !inherits(try(Revo.version,silent=TRUE),"try-error")
	if(r.open &  mac & cpu > 1)	Sys.setenv("VECLIB_MAXIMUM_THREADS" = "1")
	if(wind)	cpu <- 1
	
	cat(" Attaching data...\n")
	PHENO <- read.table(pfile, head=TRUE)
	TAXA <- colnames(PHENO)[pheno]
	PHENO <- PHENO[, pheno]
	N.Ind <- length(PHENO)
	GEBV <- PHENO
	GENO <- attach.big.matrix(paste(gfile, ".geno.desc", sep=""))
	MAP <-  try(read.table(paste(gfile, ".map", sep=""), head=FALSE), silent=TRUE)
	if((!is.null(Top.index) | !is.null(Top.perc)) & class(MAP) == "try-error"){
		stop("Please provid the Map information for all SNPs!")
	}
	if(nrow(MAP) != nrow(GENO)){
		stop("The number of SNPs in genotype and map doesn't match!")
	}
	MAP <- as.matrix(MAP)
	options(warn = -1)
	max.chr <- max(as.numeric(MAP[, 2]), na.rm=TRUE)
	if(is.infinite(max.chr))	max.chr <- 0
	map.xy.index <- which(!as.numeric(MAP[, 2]) %in% c(0 : max.chr))
	if(length(map.xy.index) != 0){
		chr.xy <- unique(MAP[map.xy.index, 2])
		for(i in 1:length(chr.xy)){
			MAP[MAP[, 2] == chr.xy[i], 2] <- max.chr + i
		}
	}
	MAP <- matrix(as.numeric(MAP), nrow(MAP))
	options(warn = 0)
	if(!is.null(kfile)){
		KIN <- read.table(kfile, head=FALSE, colClasses="numeric")
		KIN <- as.matrix(KIN)
		if(!is.numeric(KIN)){
			stop("Some none numerical values appear in KINSHIP!")
		}
		N.k <- nrow(KIN)
	}else{
		N.k <- ncol(GENO)
	}
	if(!is.null(cfile)){
		Cov <- read.table(cfile, head=FALSE, colClasses="character")
		if(sum(is.na(Cov)) != 0){
			stop("'NA' isn't allowed in 'cfile'")
		}
		Cov <- Cov[, apply(Cov, 2, function(x) length(unique(x))) != 1]
		N.c <- nrow(Cov)
		N.cc <- ncol(Cov) + 1
	}else{
		N.c <- N.Ind
		N.cc <- 1
	}
	
	cat(" Successfully attached!\n")
	
	if(!is.null(cfile)){
		Cov <- cbind(matrix(1, N.Ind), Cov)
	}else{
		Cov <- matrix(1, N.Ind)
	}
	
	cat(paste(" Trait: ", TAXA, sep=""), "\n")
	NA.index <- which(is.na(PHENO))
	NA.Ind <- length(NA.index)
	Ind.index <- c(c(1:N.Ind)[-NA.index], NA.index)
	cat(" Number of Total References:", N.Ind-NA.Ind, "\n")
	cat(" Number of Total Candidates:", NA.Ind, "\n")
	cat(" Number of Covariates:", N.cc, "\n")
	cat(" Number of Total SNPs:", nrow(GENO), "\n")
	cat(" Number of CPUs:", cpu, "\n")
	cat(" New seeds generated from:", theSeed, "\n")
	
	if((N.Ind != N.k) | (N.Ind != N.c) | (N.k != N.c)){
		stop("The numbers of individuals don't match in provided files!")
	}
	
	if(is.null(kfile)){
		cat(" Calculating marker-based Kinship...")
		KIN <- MRBLUP.Kin(GENO, type=K.method); gc()
		cat("Done!\n")
	}

	if(!is.null(prior.QTN))	Top.index <- NULL
	
	if(!is.null(Top.index) | !is.null(Top.perc)){
		cat(" Performing model: MRBLUP\n")
		cat(paste(" Cross-validation(", sample.num, "*", crv.num, ") Started...\n", sep=""))
		cross.res <- MRBLUP.CrossV(phe=PHENO, geno=GENO, prior.QTN=prior.QTN, prior.model=prior.model, vc.method=vc.method, K=KIN, map=MAP,
			max.nQTN=max.nQTN, GWAS.model=GWAS.model, theSeed=theSeed, amplify=Logx, K.method=K.method, Top.index=Top.index, bisection=bisection,
			Top.perc=Top.perc, sample.num=sample.num, crv.num=crv.num, cpu=cpu, bin.size=bin.size, bisection.loop=bisection.loop, SNP.filter=SNP.filter,
			cor.threshold=cor.threshold, judge.threshold=judge.threshold, ref.gwas=ref.gwas, GWAS.npc=GWAS.npc, CV=Cov); gc()
		cross.QTN <- cross.res$cross.QTN
		cross.model <- cross.res$cross.model
		cross.k <- cross.res$cross.k
		rm(cross.res); gc()
	}else{
		if(is.null(prior.QTN)){
			cat(" Performing model: GBLUP\n")
			cross.QTN <- NULL
			cross.model <- "K"
			cross.k <- KIN
		}else{
			cat(" Performing model: MRBLUP\n")
			cat(paste(" The provided QTNs: ", prior.QTN, "\n", sep=""))
			cat(paste(" The provided Model: ", prior.model, "\n", sep=""))
			cross.QTN <- prior.QTN
			cross.model <- prior.model
			cross.k <- KIN
		}
	}
	rm("KIN"); rm("MAP"); gc()
	cat(" Predicting...\n")
	if(is.null(cross.QTN)){
		myest <- MRBLUP.Mix(phe=PHENO, K=cross.k, vc.method=vc.method, CV=Cov)
		GEBV <- Cov %*% as.matrix(myest$beta) + myest$ebv
	}else
	{
		if(cross.model == "QTN+K"){
			if(length(cross.QTN) == 1){
				QTN.cv <- cbind(Cov, GENO[cross.QTN, ])
			}else{
				QTN.cv <- cbind(Cov, t(GENO[cross.QTN, ]))
			}
			myest <- MRBLUP.Mix(phe=PHENO, CV=QTN.cv, vc.method=vc.method, K=cross.k)
			GEBV <- myest$ebv + data.matrix(QTN.cv) %*% as.matrix(myest$beta)
		}
		if(cross.model == "QTN"){
			if(length(cross.QTN) == 1){
				QTN.cv <- as.matrix(GENO[cross.QTN, ])
			}else{
				QTN.cv <- t(GENO[cross.QTN, ])
			}
			qtn.eff <- MRBLUP.GLM(y=PHENO[-NA.index], X=Cov[-NA.index, ], qtn.matrix=QTN.cv[-NA.index, ])$QTN.eff
			GEBV <- data.matrix(cbind(Cov, QTN.cv)) %*% as.matrix(qtn.eff)
		}
	}
	#return the results
	GEBV <- as.matrix(GEBV)
	colnames(GEBV) <- TAXA

	if(file.output == TRUE){
		file.name <- paste("MRBLUP.", TAXA, ".pred", ".txt", sep="")
		write.table(GEBV, file.name, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
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
		return(paste(num, char, sep="", collapse=""))
	}
	cat(paste(" ", TAXA, " is DONE within total run time: ", times(time.cal), "\n", sep=""))
	cat(paste(c("#", rep("-", 19), "MRBLUP ACCOMPLISHED SUCCESSFULLY", rep("-", 20), "#"), collapse=""),"\n")
	return(GEBV)
}
