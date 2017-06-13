#	-----------------------------------------------------------------------------------------------
#	R package for population stratification analysis for RNA-seq (PSARseq)
#	-----------------------------------------------------------------------------------------------
#	Wrote by Gim, Jungsoo (iedenkim@gmail.com)
#	The package conducts DE analysis based on edgeR by adjusting population stratification
#	Version: 01
#	Dependencies:
#		NMF
#		edgeR


#' One-step function performing differential expression analysis adjusting population stratification
#' @param  Obj	An input object: SeqExpressionSet (EDASeq), DGEList (edgeR) or raw count in "matrix" class are available
#' @param  grps	A vector indicating the group of interest for DE test 
#' @param  n.factor	The number of factors in non-negative factorization (default = 2)
#' @details See the referenced paper.
#' @return  \item{resOut}{a list of including two variables: dgeset (DGEList) and fit (glmFit)} 
#' @export
#' @author Jungsoo Gim
#' @references Jungsoo Gim and Christoph Lange, PSARseq (2017)
#' @examples
#' cnts <- matrix(rnbinom(1000, mu=10, size=1), nrow=100)
#' grps <- factor(rep(c("a", "b"), each=5))
#' res <- PSAR(cnts, grps, n.factor=1)
#' head(DEGs(res))
PSAR <- function(Obj, grps, n.factor=2){
	if(class(Obj) == "SeqExpressionSet"){
		resOut <- PSARfromSeqExpressionSet(Obj, grps, n.factor)
	}
	else if(class(Obj) == "DGEList"){
		resOut <- PSARfromDGEList(Obj, grps, n.factor)
	}
	else if(class(Obj) == "matrix"){
		resOut <- PSARfromRawCount(Obj, grps, n.factor)
	}
	else {
		stop("object should be one of those: SeqExpressionSet, DGEList, Matrix")
	}
	return(resOut)
}

#' Extracting DEGs using LRT (edgeR) 
#' @param  PSARObj	An input object: output of PSAR() function
#' @details See the referenced paper.
#' @return  \item{DEGs}{a file including LRT result} 
#' @export
#' @author Jungsoo Gim
#' @references Jungsoo Gim and Christoph Lange, PSARseq (2017)
#' @examples
#' cnts <- matrix(rnbinom(1000, mu=10, size=1), nrow=100)
#' grps <- factor(rep(c("a", "b"), each=5))
#' res <- PSAR(cnts, grps, n.factor=1)
#' head(DEGs(res))
DEGs <- function(PSARobject){
	tmpLRT <- glmLRT(PSARobject$fit)
	tmpDEG <- topTags(tmpLRT, n = nrow(PSARobject$fit))$table
	return(tmpDEG)
}

#' Extracting covariate information including adjustment factors from non-negative matrix factorization
#' @param  PSARObj	An input object: output of PSAR() function
#' @details See the referenced paper.
#' @return  \item{dgeset}{a covariate file} 
#' @export
#' @author Jungsoo Gim
#' @references Jungsoo Gim and Christoph Lange, PSARseq (2017)
#' @examples
#' cnts <- matrix(rnbinom(1000, mu=10, size=1), nrow=100)
#' grps <- factor(rep(c("a", "b"), each=5))
#' res <- PSAR(cnts, grps, n.factor=1)
#' head(FACTORs(res))
FACTORs <- function(PSARobject){
	PSARobject$dgeset
}

PSARfromDGEList <- function(DGEset, grps, n.factor = 2){

	library(edgeR)
	library(NMF)

	if(missing(grps)){
		grps <- as.factor(DGEset$samples$group)
	}
	rawCounts <- DGEset$counts
	
	return(PSARfromRawCount(rawCounts, grps, n.factor = n.factor))

}

PSARfromSeqExpressionSet <- function(seqDataSet, grps, n.factor= 2){
	
	library(EDASeq)
	library(edgeR)
	library(NMF)
	
	if(missing(grps)){
		grps <- as.factor(pData(seqDataSet)[,1])
	}
	
	rawCounts <- counts(seqDataSet)
	return(PSARfromRawCount(rawCounts, grps, n.factor = n.factor))
}

PSARfromRawCount <- function(rawCounts, grps, n.coef=1, n.factor = 2){

	library(edgeR)
	library(NMF)
	library(RUVSeq)

	if(nrow(rawCounts) < ncol(rawCounts)){
		print("------------------------------------------------")
		print("Are you sure the input is gene by sample matrix?")
		print("------------------------------------------------")
	}
	
	if(missing(grps)){
		stop("group argument is lost: See PSARseq manual")
	}
	
	print(paste("PSARseq analysis: ", colnames(grps)[n.coef], " as a response", "\n", sep=""))
	print("Performing 1/3 step: evaluating residuals ... ")
	time1 <- Sys.time()
	designGRP <- model.matrix(~., data=grps)
	DGEset <- DGEList(counts = rawCounts, group = grps[,n.coef])
	DGEset <- calcNormFactors(DGEset, method="TMM")
	DGEset <- estimateGLMCommonDisp(DGEset, designGRP) 
	DGEset <- estimateGLMTagwiseDisp(DGEset, designGRP)
	DGEfit <- glmFit(DGEset, designGRP)
	DGEres <- residuals(DGEfit, type="deviance")
	rm(DGEset)
	
	print("Performing 2/3 step: evaluating PS factors ... ")
	resCov <- cov(DGEres)^2
	resNMF <- nmf(resCov, n.factor)
	PSARgrps <- cbind(grps, data.frame(t(coef(resNMF))))
	colnames(PSARgrps) <- c(colnames(grps), paste("H", 1:n.factor, sep="_"))
	
	print("Performing 3/3 step: Analysing DEGs ... ")
	designPSAR <- model.matrix(~., data = PSARgrps)
	DGEsetPSAR <- DGEList(counts = rawCounts, group = grps[,n.coef])
	DGEsetPSAR <- calcNormFactors(DGEsetPSAR, method="TMM")
	DGEsetPSAR <- estimateGLMCommonDisp(DGEsetPSAR, designPSAR)
	DGEsetPSAR <- estimateGLMTagwiseDisp(DGEsetPSAR, designPSAR)
	fitPSAR <- glmFit(DGEsetPSAR, designPSAR)
	LRTresult <- glmLRT(fitPSAR, coef=n.coef+1)
	time2 <- Sys.time()
	tmp.time <- time2 - time1

	print(paste("Analysis done, taking ", round(as.numeric(tmp.time),2), " ", units(tmp.time), sep=""))

	return(list(dgeset = DGEsetPSAR, fit = fitPSAR))
	
}


#	Usage example
#library(edgeR)
#library(EDASeq)
#hapmap.count <- counts(hapmap.set)
#hapmap.grps <- grps
#hapmap.DGE <- DGEList(counts = hapmap.count, group = hapmap.grps)
#rm(grps)

#set.seed(1004)
#fitCount <- PSAR(hapmap.count, hapmap.grps)
#set.seed(1004)
#fitSeqSet <- PSAR(hapmap.set)
#set.seed(1004)
#fitDGEList <- PSAR(hapmap.DGE)


#	Summarize the result
#resLRT <- glmLRT(tmpFitPSAR$fit)
#resDEG <- topTags(resLRT, n = nrow(tmpFitPSAR$fit))$table
#head(resDEG)

#	Alternative summarizing
#head(getDEGs(fitCount))
#head(getDEGs(fitSeqSet))
#head(getDEGs(fitDGEList))

