#' @docType package
#' @name MethylDeconBloodSubtypes
#' @title Deconvolute array-based methylation data.

#' Fit models.
#'
#' @family Model fitting
#' @param testVector
#' @param basis
#' @return A numeric matrix
fitNonneg.nosum1 <- function (testVector, basis){ 
	basisIndex <- c(seq(1,dim(basis)[2]))
	basisSplice <- basis
	while (1) 
	{
		afit <- lsfit (basisSplice, testVector,intercept=F)
		coef <- afit$coefficients[1:(dim(basisSplice)[2])]
		coeff <- c(rep(0,dim(basis)[2]))
		coeff[basisIndex] <- coef
		
		if (length(which(coeff < 0)) == 0) 
		{
			coeff2 <- coeff
			if (sum(coeff)>1)
			{
				coeff2 <- coeff/sum(coeff)
			}
			return(coeff2)
		} 
		
		else 
		{
			low <- which.min(coeff)
			basisIndex[low] <- 0
			basisSplice <- basis[,which(basisIndex > 0)]
			
			if ((length(which(coeff > 0)) == 1) && (length(which(coeff != 0)) == 2))
			{
				coeff[which(coeff > 0)] <- 1
				coeff[which(coeff <= 0)] <- 0
				return(coeff)
				
			}
		}
	}
}

## (may want to also hard code in this list using the default options to make code faster)
#' Get the CpG list for the model for main cell types.
#'
#' @family Generate list
#' @param length.anova An integer numeric
#' @param length.pairs An integer numeric
#' @return A character
#' @export
get.cpg.list.main <- function(length.anova = 1000, length.pairs = 300){
	## get keep list of CpGs for overall ANOVA test for all types
	ord.anova <- tests.main[,"anova.p"]
	names(ord.anova) <- row.names(tests.main)
	ord.anova <- ord.anova[order(ord.anova)]
	names.ord.anova <- names(ord.anova)
	
	## get keep list of CpGs for all pairwise comparisons
	
	## CD4-CD8
	ord.cd4.cd8 <- tests.main[,"cd4.cd8.p"]
	names(ord.cd4.cd8) <- row.names(tests.main)
	ord.cd4.cd8 <- ord.cd4.cd8[order(ord.cd4.cd8)]
	names.ord.cd4.cd8 <- names(ord.cd4.cd8)
	
	
	## CD4-CD14
	ord.cd4.cd14 <- tests.main[,"cd4.cd14.p"]
	names(ord.cd4.cd14) <- row.names(tests.main)
	ord.cd4.cd14 <- ord.cd4.cd14[order(ord.cd4.cd14)]
	names.ord.cd4.cd14 <- names(ord.cd4.cd14)
	
	
	## CD4-CD19
	ord.cd4.cd19 <- tests.main[,"cd4.cd19.p"]
	names(ord.cd4.cd19) <- row.names(tests.main)
	ord.cd4.cd19 <- ord.cd4.cd19[order(ord.cd4.cd19)]
	names.ord.cd4.cd19 <- names(ord.cd4.cd19)
	
	## CD4-gran
	ord.cd4.gran <- tests.main[,"cd4.gran.p"]
	names(ord.cd4.gran) <- row.names(tests.main)
	ord.cd4.gran <- ord.cd4.gran[order(ord.cd4.gran)]
	names.ord.cd4.gran <- names(ord.cd4.gran)
	
	
	## CD4-nk
	ord.cd4.nk <- tests.main[,"cd4.nk.p"]
	names(ord.cd4.nk) <- row.names(tests.main)
	ord.cd4.nk <- ord.cd4.nk[order(ord.cd4.nk)]
	names.ord.cd4.nk <- names(ord.cd4.nk)
	
	
	## CD8-CD14
	ord.cd8.cd14 <- tests.main[,"cd8.cd14.p"]
	names(ord.cd8.cd14) <- row.names(tests.main)
	ord.cd8.cd14 <- ord.cd8.cd14[order(ord.cd8.cd14)]
	names.ord.cd8.cd14 <- names(ord.cd8.cd14)
	
	## CD8-CD19
	ord.cd8.cd19 <- tests.main[,"cd8.cd19.p"]
	names(ord.cd8.cd19) <- row.names(tests.main)
	ord.cd8.cd19 <- ord.cd8.cd19[order(ord.cd8.cd19)]
	names.ord.cd8.cd19 <- names(ord.cd8.cd19)
	
	## CD8-gran
	ord.cd8.gran <- tests.main[,"cd8.gran.p"]
	names(ord.cd8.gran) <- row.names(tests.main)
	ord.cd8.gran <- ord.cd8.gran[order(ord.cd8.gran)]
	names.ord.cd8.gran <- names(ord.cd8.gran)
	
	## CD8-nk
	ord.cd8.nk <- tests.main[,"cd8.nk.p"]
	names(ord.cd8.nk) <- row.names(tests.main)
	ord.cd8.nk <- ord.cd8.nk[order(ord.cd8.nk)]
	names.ord.cd8.nk <- names(ord.cd8.nk)
	


	## CD14-CD19
	ord.cd14.cd19 <- tests.main[,"cd14.cd19.p"]
	names(ord.cd14.cd19) <- row.names(tests.main)
	ord.cd14.cd19 <- ord.cd14.cd19[order(ord.cd14.cd19)]
	names.ord.cd14.cd19 <- names(ord.cd14.cd19)
	
	## CD14-gran
	ord.cd14.gran <- tests.main[,"cd14.gran.p"]
	names(ord.cd14.gran) <- row.names(tests.main)
	ord.cd14.gran <- ord.cd14.gran[order(ord.cd14.gran)]
	names.ord.cd14.gran <- names(ord.cd14.gran)
	
	## CD14-nk
	ord.cd14.nk <- tests.main[,"cd14.nk.p"]
	names(ord.cd14.nk) <- row.names(tests.main)
	ord.cd14.nk <- ord.cd14.nk[order(ord.cd14.nk)]
	names.ord.cd14.nk <- names(ord.cd14.nk)
	
	
	
	## CD19-gran
	ord.cd19.gran <- tests.main[,"cd19.gran.p"]
	names(ord.cd19.gran) <- row.names(tests.main)
	ord.cd19.gran <- ord.cd19.gran[order(ord.cd19.gran)]
	names.ord.cd19.gran <- names(ord.cd19.gran)
	
	## CD19-nk
	ord.cd19.nk <- tests.main[,"cd19.nk.p"]
	names(ord.cd19.nk) <- row.names(tests.main)
	ord.cd19.nk <- ord.cd19.nk[order(ord.cd19.nk)]
	names.ord.cd19.nk <- names(ord.cd19.nk)
	
	
	## gran-nk
	ord.gran.nk <- tests.main[,"gran.nk.p"]
	names(ord.gran.nk) <- row.names(tests.main)
	ord.gran.nk <- ord.gran.nk[order(ord.gran.nk)]
	names.ord.gran.nk <- names(ord.gran.nk)
	
	
	init <- cbind(names.ord.anova,names.ord.cd4.cd8,names.ord.cd4.cd14,names.ord.cd4.cd19,names.ord.cd4.gran,names.ord.cd4.nk,names.ord.cd8.cd14,names.ord.cd8.cd19,names.ord.cd8.gran,names.ord.cd8.nk,names.ord.cd14.cd19,names.ord.cd14.gran,names.ord.cd14.nk,names.ord.cd19.gran,names.ord.cd19.nk,names.ord.gran.nk)


	keep.anova <- init[1:length.anova,"names.ord.anova"]
	
	list.rm.cd4.cd8 <- c(init[1:length.pairs,"names.ord.cd4.cd14"],init[1:length.pairs,"names.ord.cd4.cd19"],init[1:length.pairs,"names.ord.cd4.gran"],init[1:length.pairs,"names.ord.cd4.nk"],init[1:length.pairs,"names.ord.cd8.cd14"],init[1:length.pairs,"names.ord.cd8.cd19"],init[1:length.pairs,"names.ord.cd8.gran"],init[1:length.pairs,"names.ord.cd8.nk"],init[1:length.pairs,"names.ord.cd14.cd19"],init[1:length.pairs,"names.ord.cd14.gran"],init[1:length.pairs,"names.ord.cd14.nk"],init[1:length.pairs,"names.ord.cd19.gran"],init[1:length.pairs,"names.ord.cd19.nk"],init[1:length.pairs,"names.ord.gran.nk"])
	list.rm.cd4.cd8 <- unique(list.rm.cd4.cd8)
	ord.list.rm.cd4.cd8 <- setdiff(init[,"names.ord.cd4.cd8"],list.rm.cd4.cd8)
	
	list.rm.cd4.cd14 <- c(init[1:length.pairs,"names.ord.cd4.cd8"],init[1:length.pairs,"names.ord.cd4.cd19"],init[1:length.pairs,"names.ord.cd4.gran"],init[1:length.pairs,"names.ord.cd4.nk"],init[1:length.pairs,"names.ord.cd8.cd14"],init[1:length.pairs,"names.ord.cd8.cd19"],init[1:length.pairs,"names.ord.cd8.gran"],init[1:length.pairs,"names.ord.cd8.nk"],init[1:length.pairs,"names.ord.cd14.cd19"],init[1:length.pairs,"names.ord.cd14.gran"],init[1:length.pairs,"names.ord.cd14.nk"],init[1:length.pairs,"names.ord.cd19.gran"],init[1:length.pairs,"names.ord.cd19.nk"],init[1:length.pairs,"names.ord.gran.nk"])
	list.rm.cd4.cd14 <- unique(list.rm.cd4.cd14)
	ord.list.rm.cd4.cd14 <- setdiff(init[,"names.ord.cd4.cd14"],list.rm.cd4.cd14)


	list.rm.cd4.cd19 <- c(init[1:length.pairs,"names.ord.cd4.cd8"],init[1:length.pairs,"names.ord.cd4.cd14"],init[1:length.pairs,"names.ord.cd4.gran"],init[1:length.pairs,"names.ord.cd4.nk"],init[1:length.pairs,"names.ord.cd8.cd14"],init[1:length.pairs,"names.ord.cd8.cd19"],init[1:length.pairs,"names.ord.cd8.gran"],init[1:length.pairs,"names.ord.cd8.nk"],init[1:length.pairs,"names.ord.cd14.cd19"],init[1:length.pairs,"names.ord.cd14.gran"],init[1:length.pairs,"names.ord.cd14.nk"],init[1:length.pairs,"names.ord.cd19.gran"],init[1:length.pairs,"names.ord.cd19.nk"],init[1:length.pairs,"names.ord.gran.nk"])
	list.rm.cd4.cd19 <- unique(list.rm.cd4.cd19)
	ord.list.rm.cd4.cd19 <- setdiff(init[,"names.ord.cd4.cd19"],list.rm.cd4.cd19)

	list.rm.cd4.gran <- c(init[1:length.pairs,"names.ord.cd4.cd8"],init[1:length.pairs,"names.ord.cd4.cd14"],init[1:length.pairs,"names.ord.cd4.cd19"],init[1:length.pairs,"names.ord.cd4.nk"],init[1:length.pairs,"names.ord.cd8.cd14"],init[1:length.pairs,"names.ord.cd8.cd19"],init[1:length.pairs,"names.ord.cd8.gran"],init[1:length.pairs,"names.ord.cd8.nk"],init[1:length.pairs,"names.ord.cd14.cd19"],init[1:length.pairs,"names.ord.cd14.gran"],init[1:length.pairs,"names.ord.cd14.nk"],init[1:length.pairs,"names.ord.cd19.gran"],init[1:length.pairs,"names.ord.cd19.nk"],init[1:length.pairs,"names.ord.gran.nk"])
	list.rm.cd4.gran <- unique(list.rm.cd4.gran)
	ord.list.rm.cd4.gran <- setdiff(init[,"names.ord.cd4.gran"],list.rm.cd4.gran)

	list.rm.cd4.nk <- c(init[1:length.pairs,"names.ord.cd4.cd8"],init[1:length.pairs,"names.ord.cd4.cd14"],init[1:length.pairs,"names.ord.cd4.cd19"],init[1:length.pairs,"names.ord.cd4.gran"],init[1:length.pairs,"names.ord.cd8.cd14"],init[1:length.pairs,"names.ord.cd8.cd19"],init[1:length.pairs,"names.ord.cd8.gran"],init[1:length.pairs,"names.ord.cd8.nk"],init[1:length.pairs,"names.ord.cd14.cd19"],init[1:length.pairs,"names.ord.cd14.gran"],init[1:length.pairs,"names.ord.cd14.nk"],init[1:length.pairs,"names.ord.cd19.gran"],init[1:length.pairs,"names.ord.cd19.nk"],init[1:length.pairs,"names.ord.gran.nk"])
	list.rm.cd4.nk <- unique(list.rm.cd4.nk)
	ord.list.rm.cd4.nk <- setdiff(init[,"names.ord.cd4.nk"],list.rm.cd4.nk)


	list.rm.cd8.cd14 <- c(init[1:length.pairs,"names.ord.cd4.cd8"],init[1:length.pairs,"names.ord.cd4.cd14"],init[1:length.pairs,"names.ord.cd4.cd19"],init[1:length.pairs,"names.ord.cd4.gran"],init[1:length.pairs,"names.ord.cd4.nk"],init[1:length.pairs,"names.ord.cd8.cd19"],init[1:length.pairs,"names.ord.cd8.gran"],init[1:length.pairs,"names.ord.cd8.nk"],init[1:length.pairs,"names.ord.cd14.cd19"],init[1:length.pairs,"names.ord.cd14.gran"],init[1:length.pairs,"names.ord.cd14.nk"],init[1:length.pairs,"names.ord.cd19.gran"],init[1:length.pairs,"names.ord.cd19.nk"],init[1:length.pairs,"names.ord.gran.nk"])
	list.rm.cd8.cd14 <- unique(list.rm.cd8.cd14)
	ord.list.rm.cd8.cd14 <- setdiff(init[,"names.ord.cd8.cd14"],list.rm.cd8.cd14)
	
	list.rm.cd8.cd19 <- c(init[1:length.pairs,"names.ord.cd4.cd8"],init[1:length.pairs,"names.ord.cd4.cd14"],init[1:length.pairs,"names.ord.cd4.cd19"],init[1:length.pairs,"names.ord.cd4.gran"],init[1:length.pairs,"names.ord.cd4.nk"],init[1:length.pairs,"names.ord.cd8.cd14"],init[1:length.pairs,"names.ord.cd8.gran"],init[1:length.pairs,"names.ord.cd8.nk"],init[1:length.pairs,"names.ord.cd14.cd19"],init[1:length.pairs,"names.ord.cd14.gran"],init[1:length.pairs,"names.ord.cd14.nk"],init[1:length.pairs,"names.ord.cd19.gran"],init[1:length.pairs,"names.ord.cd19.nk"],init[1:length.pairs,"names.ord.gran.nk"])
	list.rm.cd8.cd19 <- unique(list.rm.cd8.cd19)
	ord.list.rm.cd8.cd19 <- setdiff(init[,"names.ord.cd8.cd19"],list.rm.cd8.cd19)
	
	list.rm.cd8.gran <- c(init[1:length.pairs,"names.ord.cd4.cd8"],init[1:length.pairs,"names.ord.cd4.cd14"],init[1:length.pairs,"names.ord.cd4.cd19"],init[1:length.pairs,"names.ord.cd4.gran"],init[1:length.pairs,"names.ord.cd4.nk"],init[1:length.pairs,"names.ord.cd8.cd14"],init[1:length.pairs,"names.ord.cd8.cd19"],init[1:length.pairs,"names.ord.cd8.nk"],init[1:length.pairs,"names.ord.cd14.cd19"],init[1:length.pairs,"names.ord.cd14.gran"],init[1:length.pairs,"names.ord.cd14.nk"],init[1:length.pairs,"names.ord.cd19.gran"],init[1:length.pairs,"names.ord.cd19.nk"],init[1:length.pairs,"names.ord.gran.nk"])
	list.rm.cd8.gran <- unique(list.rm.cd8.gran)
	ord.list.rm.cd8.gran <- setdiff(init[,"names.ord.cd8.gran"],list.rm.cd8.gran)
	
	
	list.rm.cd8.nk <- c(init[1:length.pairs,"names.ord.cd4.cd8"],init[1:length.pairs,"names.ord.cd4.cd14"],init[1:length.pairs,"names.ord.cd4.cd19"],init[1:length.pairs,"names.ord.cd4.gran"],init[1:length.pairs,"names.ord.cd4.nk"],init[1:length.pairs,"names.ord.cd8.cd14"],init[1:length.pairs,"names.ord.cd8.cd19"],init[1:length.pairs,"names.ord.cd8.gran"],init[1:length.pairs,"names.ord.cd14.cd19"],init[1:length.pairs,"names.ord.cd14.gran"],init[1:length.pairs,"names.ord.cd14.nk"],init[1:length.pairs,"names.ord.cd19.gran"],init[1:length.pairs,"names.ord.cd19.nk"],init[1:length.pairs,"names.ord.gran.nk"])
	list.rm.cd8.nk <- unique(list.rm.cd8.nk)
	ord.list.rm.cd8.nk <- setdiff(init[,"names.ord.cd8.nk"],list.rm.cd8.nk)

	list.rm.cd14.cd19 <- c(init[1:length.pairs,"names.ord.cd4.cd8"],init[1:length.pairs,"names.ord.cd4.cd14"],init[1:length.pairs,"names.ord.cd4.cd19"],init[1:length.pairs,"names.ord.cd4.gran"],init[1:length.pairs,"names.ord.cd4.nk"],init[1:length.pairs,"names.ord.cd8.cd14"],init[1:length.pairs,"names.ord.cd8.cd19"],init[1:length.pairs,"names.ord.cd8.gran"],init[1:length.pairs,"names.ord.cd8.nk"],init[1:length.pairs,"names.ord.cd14.gran"],init[1:length.pairs,"names.ord.cd14.nk"],init[1:length.pairs,"names.ord.cd19.gran"],init[1:length.pairs,"names.ord.cd19.nk"],init[1:length.pairs,"names.ord.gran.nk"])
	list.rm.cd14.cd19 <- unique(list.rm.cd14.cd19)
	ord.list.rm.cd14.cd19 <- setdiff(init[,"names.ord.cd14.cd19"],list.rm.cd14.cd19)
	
	
	list.rm.cd14.gran <- c(init[1:length.pairs,"names.ord.cd4.cd8"],init[1:length.pairs,"names.ord.cd4.cd14"],init[1:length.pairs,"names.ord.cd4.cd19"],init[1:length.pairs,"names.ord.cd4.gran"],init[1:length.pairs,"names.ord.cd4.nk"],init[1:length.pairs,"names.ord.cd8.cd14"],init[1:length.pairs,"names.ord.cd8.cd19"],init[1:length.pairs,"names.ord.cd8.gran"],init[1:length.pairs,"names.ord.cd8.nk"],init[1:length.pairs,"names.ord.cd14.cd19"],init[1:length.pairs,"names.ord.cd14.nk"],init[1:length.pairs,"names.ord.cd19.gran"],init[1:length.pairs,"names.ord.cd19.nk"],init[1:length.pairs,"names.ord.gran.nk"])
	list.rm.cd14.gran <- unique(list.rm.cd14.gran)
	ord.list.rm.cd14.gran <- setdiff(init[,"names.ord.cd14.gran"],list.rm.cd14.gran)
	
	
	list.rm.cd14.nk <- c(init[1:length.pairs,"names.ord.cd4.cd8"],init[1:length.pairs,"names.ord.cd4.cd14"],init[1:length.pairs,"names.ord.cd4.cd19"],init[1:length.pairs,"names.ord.cd4.gran"],init[1:length.pairs,"names.ord.cd4.nk"],init[1:length.pairs,"names.ord.cd8.cd14"],init[1:length.pairs,"names.ord.cd8.cd19"],init[1:length.pairs,"names.ord.cd8.gran"],init[1:length.pairs,"names.ord.cd8.nk"],init[1:length.pairs,"names.ord.cd14.cd19"],init[1:length.pairs,"names.ord.cd14.gran"],init[1:length.pairs,"names.ord.cd19.gran"],init[1:length.pairs,"names.ord.cd19.nk"],init[1:length.pairs,"names.ord.gran.nk"])
	list.rm.cd14.nk <- unique(list.rm.cd14.nk)
	ord.list.rm.cd14.nk <- setdiff(init[,"names.ord.cd14.nk"],list.rm.cd14.nk)


	list.rm.cd19.gran <- c(init[1:length.pairs,"names.ord.cd4.cd8"],init[1:length.pairs,"names.ord.cd4.cd14"],init[1:length.pairs,"names.ord.cd4.cd19"],init[1:length.pairs,"names.ord.cd4.gran"],init[1:length.pairs,"names.ord.cd4.nk"],init[1:length.pairs,"names.ord.cd8.cd14"],init[1:length.pairs,"names.ord.cd8.cd19"],init[1:length.pairs,"names.ord.cd8.gran"],init[1:length.pairs,"names.ord.cd8.nk"],init[1:length.pairs,"names.ord.cd14.cd19"],init[1:length.pairs,"names.ord.cd14.gran"],init[1:length.pairs,"names.ord.cd14.nk"],init[1:length.pairs,"names.ord.cd19.nk"],init[1:length.pairs,"names.ord.gran.nk"])
	list.rm.cd19.gran <- unique(list.rm.cd19.gran)
	ord.list.rm.cd19.gran <- setdiff(init[,"names.ord.cd19.gran"],list.rm.cd19.gran)



	list.rm.cd19.nk <- c(init[1:length.pairs,"names.ord.cd4.cd8"],init[1:length.pairs,"names.ord.cd4.cd14"],init[1:length.pairs,"names.ord.cd4.cd19"],init[1:length.pairs,"names.ord.cd4.gran"],init[1:length.pairs,"names.ord.cd4.nk"],init[1:length.pairs,"names.ord.cd8.cd14"],init[1:length.pairs,"names.ord.cd8.cd19"],init[1:length.pairs,"names.ord.cd8.gran"],init[1:length.pairs,"names.ord.cd8.nk"],init[1:length.pairs,"names.ord.cd14.cd19"],init[1:length.pairs,"names.ord.cd14.gran"],init[1:length.pairs,"names.ord.cd14.nk"],init[1:length.pairs,"names.ord.cd19.gran"],init[1:length.pairs,"names.ord.gran.nk"])
	list.rm.cd19.nk <- unique(list.rm.cd19.nk)
	ord.list.rm.cd19.nk <- setdiff(init[,"names.ord.cd19.nk"],list.rm.cd19.nk)


	list.rm.gran.nk <- c(init[1:length.pairs,"names.ord.cd4.cd8"],init[1:length.pairs,"names.ord.cd4.cd14"],init[1:length.pairs,"names.ord.cd4.cd19"],init[1:length.pairs,"names.ord.cd4.gran"],init[1:length.pairs,"names.ord.cd4.nk"],init[1:length.pairs,"names.ord.cd8.cd14"],init[1:length.pairs,"names.ord.cd8.cd19"],init[1:length.pairs,"names.ord.cd8.gran"],init[1:length.pairs,"names.ord.cd8.nk"],init[1:length.pairs,"names.ord.cd14.cd19"],init[1:length.pairs,"names.ord.cd14.gran"],init[1:length.pairs,"names.ord.cd14.nk"],init[1:length.pairs,"names.ord.cd19.gran"],init[1:length.pairs,"names.ord.cd19.nk"])
	list.rm.gran.nk <- unique(list.rm.gran.nk)
	ord.list.rm.gran.nk <- setdiff(init[,"names.ord.gran.nk"],list.rm.gran.nk)


	keep.cd4.cd8 <- ord.list.rm.cd4.cd8[1:length.pairs]
	keep.cd4.cd14 <- ord.list.rm.cd4.cd14[1:length.pairs]
	keep.cd4.cd19 <- ord.list.rm.cd4.cd19[1:length.pairs]
	keep.cd4.gran <- ord.list.rm.cd4.gran[1:length.pairs]
	keep.cd4.nk <- ord.list.rm.cd4.nk[1:length.pairs]
	
	keep.cd8.cd14 <- ord.list.rm.cd8.cd14[1:length.pairs]
	keep.cd8.cd19 <- ord.list.rm.cd8.cd19[1:length.pairs]
	keep.cd8.gran <- ord.list.rm.cd8.gran[1:length.pairs]
	keep.cd8.nk <- ord.list.rm.cd8.nk[1:length.pairs]
	
	keep.cd14.cd19 <- ord.list.rm.cd14.cd19[1:length.pairs]
	keep.cd14.gran <- ord.list.rm.cd14.gran[1:length.pairs]
	keep.cd14.nk <- ord.list.rm.cd14.nk[1:length.pairs]
	
	keep.cd19.gran <- ord.list.rm.cd19.gran[1:length.pairs]
	keep.cd19.nk <- ord.list.rm.cd19.nk[1:length.pairs]
	
	keep.gran.nk <- ord.list.rm.gran.nk[1:length.pairs]	

	
	## merge CpG sets to get final cpg list
	cpg.list <- c(keep.anova,keep.cd4.cd8,keep.cd4.cd14,keep.cd4.cd19,keep.cd4.gran,keep.cd4.nk,keep.cd8.cd14,keep.cd8.cd19,keep.cd8.gran,keep.cd8.nk,keep.cd14.cd19,keep.cd14.gran,keep.cd14.nk,keep.cd19.gran,keep.cd19.nk,keep.gran.nk)
	
	cpg.list <- unique(cpg.list)
	
	return(cpg.list)
}

#' Fit the model for the 6 main cell types
#'
#' @family Model fitting
#' @param data A data frame
#' @param cpg.list A character
#' @return A character
#' @export
fit.main.nosum1 <- function(data, cpg.list = cpg.list.main.stage1){

	cpg.list <- intersect(cpg.list,row.names(data))
		
	fit.sub3 <- apply(data[match(cpg.list,row.names(data)),],2,fitNonneg.nosum1,med.main[match(cpg.list,row.names(med.main)),])

	fit.sub3 <- t(fit.sub3)
	

	colnames(fit.sub3) <- c("CD4","CD8","CD14","CD19","Gran","NK")

	return(fit.sub3)
}

## Functions to partition the whole blood beta in order to fit the two-stage model (refine CD4/CD8 estimates) 
## 	and also to fit T and B cell subtype models
#' Partition whole blood beta for two-stage model fit inner
#'
#' @family Data partitioning
#' @param x
#' @param fit
#' @param keep.cols
#' @return A numeric matrix
get.beta.sub.inner <- function(x, fit, keep.cols){
	x <- as.numeric(x)
	keep.len <- length(keep.cols)
	sub.cols <- setdiff(seq(from=1,to=dim(fit)[2]),keep.cols)
	sub.len <- length(sub.cols)
	beta <- x[(sub.len+1):length(x)]
	blood <- x[1:sub.len]
	
	temp.sum <- 0
	for (i in (1:sub.len))
	{
		temp <- fit[,sub.cols[i]]*as.numeric(blood[i])
		temp.sum <- temp.sum+temp
	}
	temp2 <- beta-temp.sum
	return(temp2)
}

#' Partition whole blood beta for two-stage model fit
#' 
#' @family Data partitioning
#' @param data A data frame
#' @param blood
#' @param fit
#' @param keep.cols
#' @return A numeric matrix
#' @export
get.beta.sub <- function(data, blood = med.main, fit, keep.cols){
	per.sum <- apply(fit,1,function(x){sum(x[keep.cols])})
	#blood.sub <- blood[,keep.cols]
	sub.cols <- setdiff(seq(from=1,to=dim(fit)[2]),keep.cols)
	keep.cpgs <- intersect(row.names(data),row.names(blood))
	data.all <- cbind(blood[keep.cpgs,sub.cols],data[keep.cpgs,])
	beta.sub <- t(apply(data.all,1,get.beta.sub.inner,fit=fit,keep.cols=keep.cols))
	beta.sub2 <- beta.sub/per.sum
	return(beta.sub)
}

## Stage 2 models to improve CD4/CD8 main cell type estimates -------------
## Stage 2b (CD4 vs. CD8)

# (should hard code in CpG list also)
#' Get CpG list for stage 2
#' 
#' @family Generate list
#' @param length.pairs An integer numeric
#' @return A Character
#' @export
get.cpg.list.stage2b <- function(length.pairs = 100){
	## get keep list of CpGs for all pairwise comparisons
	
	## CD4-CD8
	ord.cd4.cd8 <- tests.main[,"cd4.cd8.p"]
	names(ord.cd4.cd8) <- row.names(tests.main)
	ord.cd4.cd8 <- ord.cd4.cd8[order(ord.cd4.cd8)]
	names.ord.cd4.cd8 <- names(ord.cd4.cd8)
	
	
	## CD4-CD14
	ord.cd4.cd14 <- tests.main[,"cd4.cd14.p"]
	names(ord.cd4.cd14) <- row.names(tests.main)
	ord.cd4.cd14 <- ord.cd4.cd14[order(ord.cd4.cd14)]
	names.ord.cd4.cd14 <- names(ord.cd4.cd14)
	
	
	## CD4-CD19
	ord.cd4.cd19 <- tests.main[,"cd4.cd19.p"]
	names(ord.cd4.cd19) <- row.names(tests.main)
	ord.cd4.cd19 <- ord.cd4.cd19[order(ord.cd4.cd19)]
	names.ord.cd4.cd19 <- names(ord.cd4.cd19)
	
	## CD4-gran
	ord.cd4.gran <- tests.main[,"cd4.gran.p"]
	names(ord.cd4.gran) <- row.names(tests.main)
	ord.cd4.gran <- ord.cd4.gran[order(ord.cd4.gran)]
	names.ord.cd4.gran <- names(ord.cd4.gran)
	
	
	## CD4-nk
	ord.cd4.nk <- tests.main[,"cd4.nk.p"]
	names(ord.cd4.nk) <- row.names(tests.main)
	ord.cd4.nk <- ord.cd4.nk[order(ord.cd4.nk)]
	names.ord.cd4.nk <- names(ord.cd4.nk)
	
	
	## CD8-CD14
	ord.cd8.cd14 <- tests.main[,"cd8.cd14.p"]
	names(ord.cd8.cd14) <- row.names(tests.main)
	ord.cd8.cd14 <- ord.cd8.cd14[order(ord.cd8.cd14)]
	names.ord.cd8.cd14 <- names(ord.cd8.cd14)
	
	## CD8-CD19
	ord.cd8.cd19 <- tests.main[,"cd8.cd19.p"]
	names(ord.cd8.cd19) <- row.names(tests.main)
	ord.cd8.cd19 <- ord.cd8.cd19[order(ord.cd8.cd19)]
	names.ord.cd8.cd19 <- names(ord.cd8.cd19)
	
	## CD8-gran
	ord.cd8.gran <- tests.main[,"cd8.gran.p"]
	names(ord.cd8.gran) <- row.names(tests.main)
	ord.cd8.gran <- ord.cd8.gran[order(ord.cd8.gran)]
	names.ord.cd8.gran <- names(ord.cd8.gran)
	
	## CD8-nk
	ord.cd8.nk <- tests.main[,"cd8.nk.p"]
	names(ord.cd8.nk) <- row.names(tests.main)
	ord.cd8.nk <- ord.cd8.nk[order(ord.cd8.nk)]
	names.ord.cd8.nk <- names(ord.cd8.nk)
	


	## CD14-CD19
	ord.cd14.cd19 <- tests.main[,"cd14.cd19.p"]
	names(ord.cd14.cd19) <- row.names(tests.main)
	ord.cd14.cd19 <- ord.cd14.cd19[order(ord.cd14.cd19)]
	names.ord.cd14.cd19 <- names(ord.cd14.cd19)
	
	## CD14-gran
	ord.cd14.gran <- tests.main[,"cd14.gran.p"]
	names(ord.cd14.gran) <- row.names(tests.main)
	ord.cd14.gran <- ord.cd14.gran[order(ord.cd14.gran)]
	names.ord.cd14.gran <- names(ord.cd14.gran)
	
	## CD14-nk
	ord.cd14.nk <- tests.main[,"cd14.nk.p"]
	names(ord.cd14.nk) <- row.names(tests.main)
	ord.cd14.nk <- ord.cd14.nk[order(ord.cd14.nk)]
	names.ord.cd14.nk <- names(ord.cd14.nk)
	
	
	
	## CD19-gran
	ord.cd19.gran <- tests.main[,"cd19.gran.p"]
	names(ord.cd19.gran) <- row.names(tests.main)
	ord.cd19.gran <- ord.cd19.gran[order(ord.cd19.gran)]
	names.ord.cd19.gran <- names(ord.cd19.gran)
	
	## CD19-nk
	ord.cd19.nk <- tests.main[,"cd19.nk.p"]
	names(ord.cd19.nk) <- row.names(tests.main)
	ord.cd19.nk <- ord.cd19.nk[order(ord.cd19.nk)]
	names.ord.cd19.nk <- names(ord.cd19.nk)
	
	
	## gran-nk
	ord.gran.nk <- tests.main[,"gran.nk.p"]
	names(ord.gran.nk) <- row.names(tests.main)
	ord.gran.nk <- ord.gran.nk[order(ord.gran.nk)]
	names.ord.gran.nk <- names(ord.gran.nk)
	
	
	init <- cbind(names.ord.cd4.cd8,names.ord.cd4.cd14,names.ord.cd4.cd19,names.ord.cd4.gran,names.ord.cd4.nk,names.ord.cd8.cd14,names.ord.cd8.cd19,names.ord.cd8.gran,names.ord.cd8.nk,names.ord.cd14.cd19,names.ord.cd14.gran,names.ord.cd14.nk,names.ord.cd19.gran,names.ord.cd19.nk,names.ord.gran.nk)

	
	list.rm.cd4.cd8 <- c(init[1:length.pairs,"names.ord.cd4.cd14"],init[1:length.pairs,"names.ord.cd4.cd19"],init[1:length.pairs,"names.ord.cd4.gran"],init[1:length.pairs,"names.ord.cd4.nk"],init[1:length.pairs,"names.ord.cd8.cd14"],init[1:length.pairs,"names.ord.cd8.cd19"],init[1:length.pairs,"names.ord.cd8.gran"],init[1:length.pairs,"names.ord.cd8.nk"],init[1:length.pairs,"names.ord.cd14.cd19"],init[1:length.pairs,"names.ord.cd14.gran"],init[1:length.pairs,"names.ord.cd14.nk"],init[1:length.pairs,"names.ord.cd19.gran"],init[1:length.pairs,"names.ord.cd19.nk"],init[1:length.pairs,"names.ord.gran.nk"])
	list.rm.cd4.cd8 <- unique(list.rm.cd4.cd8)
	ord.list.rm.cd4.cd8 <- setdiff(init[,"names.ord.cd4.cd8"],list.rm.cd4.cd8)	
	
	
	
	
	keep.cd4.cd8 <- ord.list.rm.cd4.cd8[1:length.pairs]

	
	## merge CpG sets to get final cpg list
	cpg.list <- keep.cd4.cd8
	
	cpg.list <- unique(cpg.list)
	
	return(cpg.list)
}

#' Fit models for stage 2
#' 
#' @family Model fitting
#' @param data A data frame
#' @param fit.main A numeric matrix
#' @param cpg.list A character
#' @return A numeric matrix
fit.stage2b <- function(data, fit.main, cpg.list = cpg.list.2b){

	t.per <- fit.main[,1]+fit.main[,2]
	cpg.list <- intersect(cpg.list,row.names(data))
	
	med.t <- med.main[,c(1,2)]
		
	fit.sub3 <- apply(data[match(cpg.list,row.names(data)),],2,fitNonneg.nosum1,med.t[match(cpg.list,row.names(med.t)),])

	fit.sub3 <- t(fit.sub3)
	
	sum.fit <- apply(fit.sub3,1,sum)
	
	ratio <- sum.fit/t.per
	
	fit.sub3 <- fit.sub3/ratio	
	
	fit.new <- cbind(fit.sub3[,1],fit.sub3[,2],fit.main[,3],fit.main[,4],fit.main[,5],fit.main[,6])
	colnames(fit.new) <- c("CD4","CD8","CD14","CD19","Gran","NK")

	return(fit.new)
}


### functions for CD4 subtypes -----------------------------------------------

## get CpG list for CD4 subtypes (should hard code in CpG list also)
#' Get CpG list for CD4 subtypes
#' @family Generate list
#' @param length.anova An integer numeric
#' @param length.pairs An integer numeric
#' @return A character
get.cpg.list.cd4.sub <- function(length.anova = 2100, length.pairs = 600){
	## get keep list of CpGs for overall ANOVA test for all types
	ord.anova <- tests.cd4.sub[,"cd4.anova.p"]
	names(ord.anova) <- row.names(tests.cd4.sub)
	ord.anova <- ord.anova[order(ord.anova)]
	names.ord.anova <- names(ord.anova)
	
	## get keep list of CpGs for all pairwise comparisons
	
	## Mem-Naive
	ord.mem.naive <- tests.cd4.sub[,"cd4.mem.naive.p"]
	names(ord.mem.naive) <- row.names(tests.cd4.sub)
	ord.mem.naive <- ord.mem.naive[order(ord.mem.naive)]
	names.ord.mem.naive <- names(ord.mem.naive)
	
	
	## Mem-Reg
	ord.mem.reg <- tests.cd4.sub[,"cd4.mem.reg.p"]
	names(ord.mem.reg) <- row.names(tests.cd4.sub)
	ord.mem.reg <- ord.mem.reg[order(ord.mem.reg)]
	names.ord.mem.reg <- names(ord.mem.reg)
	
	
	## Naive-Reg
	ord.naive.reg <- tests.cd4.sub[,"cd4.naive.reg.p"]
	names(ord.naive.reg) <- row.names(tests.cd4.sub)
	ord.naive.reg <- ord.naive.reg[order(ord.naive.reg)]
	names.ord.naive.reg <- names(ord.naive.reg)


	keep.anova <- names.ord.anova[1:length.anova]
	
	list.rm.mem.naive <- c(names.ord.mem.reg[1:length.pairs],names.ord.naive.reg[1:length.pairs])
	list.rm.mem.naive <- unique(list.rm.mem.naive)
	ord.list.rm.mem.naive <- setdiff(names.ord.mem.naive,list.rm.mem.naive)
	
	list.rm.mem.reg <- c(names.ord.mem.naive[1:length.pairs],names.ord.naive.reg[1:length.pairs])
	list.rm.mem.reg <- unique(list.rm.mem.reg)
	ord.list.rm.mem.reg <- setdiff(names.ord.mem.reg,list.rm.mem.reg)
	
	list.rm.naive.reg <- c(names.ord.mem.naive[1:length.pairs],names.ord.mem.reg[1:length.pairs])
	list.rm.naive.reg <- unique(list.rm.naive.reg)
	ord.list.rm.naive.reg <- setdiff(names.ord.naive.reg,list.rm.naive.reg)
	
	
	
	keep.mem.naive <- ord.list.rm.mem.naive[1:length.pairs]
	keep.mem.reg <- ord.list.rm.mem.reg[1:length.pairs]
	keep.naive.reg <- ord.list.rm.naive.reg[1:length.pairs]
		
	## merge CpG sets to get final cpg list
	cpg.list <- c(keep.anova,keep.mem.naive,keep.mem.reg,keep.naive.reg)
	
	cpg.list <- unique(cpg.list)
	
	return(cpg.list)
}

#(need to partition data into CD4 beta value first)
#' Fit model for CD4 subtypes on whole blood data
#' 
#' @family Generate list
#' @param data A data frame
#' @param fit.main A numeric matrix
#' @param cpg.list A character
#' @return A numeric matrix
#' @export
fit.cd4.sub <- function(data, fit.main, cpg.list = cpg.list.cd4) {

	cd4.per <- fit.main[,1]
	cpg.list <- intersect(cpg.list,row.names(data))
		
	fit.sub3 <- apply(data[match(cpg.list,row.names(data)),],2,fitNonneg.nosum1,med.cd4.sub[match(cpg.list,row.names(med.cd4.sub)),])

	fit.sub3 <- t(fit.sub3)
	
	sum.fit <- apply(fit.sub3,1,sum)
	
	ratio <- sum.fit/cd4.per
	
	fit.sub3 <- fit.sub3/ratio	
	colnames(fit.sub3) <- c("CD4-Naive","CD4-Mem","CD4-Reg")

	return(fit.sub3)
}

#' Fit model for CD4 subtypes on M450 data from sorted CD4 T cells
#' 
#' @family Model fitting
#' @param data A data frame
#' @param cpg.list A character
#' @return A numeric matrix
fit.cd4.sorted <- function(data, cpg.list = cpg.list.cd4){
	#cd4.per <- fit.main[,1]
	cpg.list <- intersect(cpg.list,row.names(data))
		
	fit.sub3 <- apply(data[match(cpg.list,row.names(data)),],2,fitNonneg.nosum1,med.cd4.sub[match(cpg.list,row.names(med.cd4.sub)),])

	fit.sub3 <- t(fit.sub3)
	
	#sum.fit <- apply(fit.sub3,1,sum)
	
	#ratio <- sum.fit/cd4.per
	
	#fit.sub3 <- fit.sub3/ratio	
	colnames(fit.sub3) <- c("CD4-Naive","CD4-Mem","CD4-Reg")

	return(fit.sub3)
}


## functions to fit model for CD8 subtypes ------------------------------------

#  (should hard code in this CpG list also)
#' Get CpG list for CD8 subtypes
#' 
#' @family Generate list
#' @param length.list An integer numeric
#' @return A character
#' @export
get.cpg.list.cd8.sub <- function(length.list = 300){
	## get keep list of CpGs for overall ANOVA test for all types
	ord.list <- tests.cd8.sub[,"cd8.mem.naive.p"]
	names(ord.list) <- row.names(tests.cd8.sub)
	ord.list <- ord.list[order(ord.list)]
	names.ord.list <- names(ord.list)
	
		## merge CpG sets to get final cpg list
	cpg.list <- names.ord.list[1:length.list]
	
	cpg.list <- unique(cpg.list)
	
	return(cpg.list)
}

# (need to parition data into CD8 beta value first)
#' Fit model for CD8 subtypes for whole blood data
#' 
#' @family Model fitting
#' @param data A data frame
#' @param fit.main A numeric matrix
#' @param cpg.list A character
#' @return A numeric matrix
#' @export
fit.cd8.sub <- function(data, fit.main, cpg.list = cpg.list.cd8){

	cd8.per <- fit.main[,2]
	cpg.list <- intersect(cpg.list,row.names(data))
		
	fit.sub3 <- apply(data[match(cpg.list,row.names(data)),],2,fitNonneg.nosum1,med.cd8.sub[match(cpg.list,row.names(med.cd8.sub)),])

	fit.sub3 <- t(fit.sub3)
	
	sum.fit <- apply(fit.sub3,1,sum)
	
	ratio <- sum.fit/cd8.per
	
	fit.sub3 <- fit.sub3/ratio	
	colnames(fit.sub3) <- c("CD8-Naive","CD8-Mem")

	return(fit.sub3)
}

#' Fit model for CD8 subtypes using M450 data from sorted CD8 T cells
#' 
#' @family Model fitting
#' @param data A data frame
#' @param cpg.list A character
#' @return A numeric matrix
fit.cd8.sub.sorted <- function(data, cpg.list = cpg.list.cd8){

	#cd8.per <- fit.main[,2]
	cpg.list <- intersect(cpg.list,row.names(data))
		
	fit.sub3 <- apply(data[match(cpg.list,row.names(data)),],2,fitNonneg.nosum1,med.cd8.sub[match(cpg.list,row.names(med.cd8.sub)),])

	fit.sub3 <- t(fit.sub3)
	
	#sum.fit <- apply(fit.sub3,1,sum)
	
	#ratio <- sum.fit/cd8.per
	
	#fit.sub3 <- fit.sub3/ratio	
	colnames(fit.sub3) <- c("CD8-Naive","CD8-Mem")

	return(fit.sub3)
}

## models for CD19 B cell subtypes -----------------------------------------------

#(should hard code in CpG list also)
#' Function to get CpG list for CD19 B cell subtypes 
#' 
#' @family Model fitting
#' @param length.pairs An integer numeric
#' @return A character
#' @export
get.cpg.list.cd19.sub <- function(length.pairs = 300){
	## get keep list of CpGs for overall ANOVA test for all types
	ord.anova <- tests.cd19.sub[,"cd19.naive.mem.p"]
	names(ord.anova) <- row.names(tests.cd19.sub)
	ord.anova <- ord.anova[order(ord.anova)]
	names.ord.anova <- names(ord.anova)
	
	
	keep.anova <- names.ord.anova[1:length.pairs]
	
	cpg.list <- keep.anova
	
	cpg.list <- unique(cpg.list)
	
	return(cpg.list)
}

# (need to partition data into B cell beta score first)
#!! med.cd19.sub not available !!
#' Fit model for CD19 B cell subtypes for whole blood  
#' 
#' @family Model fitting
#' @param data A data frame
#' @param fit.main A numeric matrix
#' @param cpg.list A character
#' @return A numeric matrix
#' @export
fit.cd19.sub <- function(data, fit.main, cpg.list = cpg.list.cd19){

	cd19.per <- fit.main[,4]
	
	cpg.list <- intersect(cpg.list,row.names(data))
		
	fit.sub3 <- apply(data[match(cpg.list,row.names(data)),],2,fitNonneg.nosum1,med.cd19.sub.2class[match(cpg.list,row.names(med.cd19.sub.2class)),])

	fit.sub3 <- t(fit.sub3)
	
	sum.fit <- apply(fit.sub3,1,sum)
	
	ratio <- sum.fit/cd19.per
	
	fit.sub3 <- fit.sub3/ratio
	
	colnames(fit.sub3) <- c("CD19-Naive","CD19-Mem")

	return(fit.sub3)
}


#' Fit model for CD19 subtypes using M450 data from sorted B cells
#' 
#' @family Model fitting
#' @param data A data frame
#' @param cpg.list A character
#' @return A numeric matrix
#' @export
fit.cd19.sorted <- function(data, cpg.list = cpg.list.cd19){

	cpg.list <- intersect(cpg.list,row.names(data))
		
	fit.sub3 <- apply(data[match(cpg.list,row.names(data)),],2,fitNonneg.nosum1,med.cd19.sub.2class[match(cpg.list,row.names(med.cd19.sub.2class)),])

	fit.sub3 <- t(fit.sub3)
	
	colnames(fit.sub3) <- c("CD19-Naive","CD19-Mem")

	return(fit.sub3)
}

#' Get stage 1 & 2 estimates and T and B cell subtypes
#' 
#' @family Subtype estimation
#' @param data A data frame
#' @return A numeric matrix
#' @export 
est.all.wb <- function(data){
	## stage 1
	fit.stage1 <- fit.main.nosum1(data)
	
	## stage 2
	# t.data <- get.beta.sub(data[cpg.list.2b,],med.main,fit.stage1,c(1,2))
	# fit.stage2 <- fit.stage2b(t.data,fit.stage1,cpg.list.2b)
	
	## cd4 subtypes
	# data.cd4.partition <- get.beta.sub(data[cpg.list.cd4,],med.main,fit.stage2,1)
	# fit.data.cd4.sub <- fit.cd4.sub(data.cd4.partition,fit.stage2,cpg.list.cd4)

	# ## cd8 subtypes
	# data.cd8.partition <- get.beta.sub(data[cpg.list.cd8,],med.main,fit.stage2,2)
	# fit.data.cd8.sub <- fit.cd8.sub(data.cd8.partition,fit.stage2,cpg.list.cd8)

	# ## cd19 subtypes
	# data.cd19.partition <- get.beta.sub(data[cpg.list.cd19,],med.main,fit.stage2,4)
	# fit.data.cd19.sub <- fit.cd19.sub(data.cd19.partition,fit.stage2,cpg.list.cd19)

	# ## merge results together
	# res <- cbind(fit.stage2,fit.data.cd4.sub,fit.data.cd8.sub,fit.data.cd19.sub)
	
	# return(res)
}
