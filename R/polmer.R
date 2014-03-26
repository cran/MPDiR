polmer <- function(formula, data, lnk = "logit", which.lme4 = "lme4.0",  ...){
	if(require(which.lme4, character.only = TRUE)) invisible() else stop("need lme4 package installed!\n")
# Get fixed and random terms
   fxForm <- lme4.0:::nobars(formula)
   fxForm[[2]] <- NULL
   rTerms <- lme4.0:::findbars(formula)
   ranTerms <- sapply(rTerms, deparse)
   fixTerms <- labels(terms(fxForm))
	ranNames <- as.vector(sapply(ranTerms, function(x)
		strsplit(x, "\\| ")[[1]][2]
		))
	LeftRanNames <- as.vector(sapply(ranTerms, function(x)
		strsplit(x, " | ", fixed = TRUE)[[1]][1]
		))
	NoRanEff <- !(length(ranTerms) > 0)

		
# Make fixed-effect model matrix
	Resp <- ordered(data[[as.character(formula[[2]])]])
	Cuts <- as.numeric(levels(Resp)[-length(levels(Resp))])
	cumRat <- as.vector(sapply(Cuts, 
		function(x) Resp <= x))
	fRat <- gl(length(Cuts), nrow(data), 
		nrow(data) * length(Cuts))
	X <- model.matrix(~ fRat - 1)
	labs <- sapply(Cuts, function(x)
		paste(x, "|", x+1, sep = "")
		)
	colnames(X) <- labs
	
	fX <- -model.matrix(fxForm, data)[, -1]
	fX.names <- if (inherits(fX, "matrix")) colnames(fX) else paste(fixTerms, levels(data[[fixTerms]])[2], sep = "")
	fX <- kronecker(matrix(rep(1, length(Cuts)), ncol = 1),
		 fX)
	X <- cbind(X, fX)
	colnames(X)[-seq(length(Cuts))] <- fX.names
	p.df <- data.frame(cumRat = cumRat, X = X)
	names(p.df) <- c("cumRat", colnames(X))

# Make random-effect model vectors
	frm <- if(!NoRanEff){
		tmp <- sapply(seq_len(length(ranNames)), 
			function(x)
				rep(data[[ranNames[x]]], length(Cuts)))

		for(ix in seq_len(ncol(tmp))) 
			assign(ranNames[ix], tmp[, ix])
		rxForm <- paste(paste("(", ranTerms, ")", 
			sep = "", collapse = " + "), " - 1")
	as.formula(paste("cumRat ~ .  + ", rxForm))
	} else as.formula("cumRat ~ . - 1")
	for(ix in LeftRanNames)
	  if(ix != "1") 
	  	assign(ix, rep(data[[ix]], each = length(Cuts)))
	  	
	CLL <- if (NoRanEff)
	  substitute(
	  glm(FRM, data = p.df, 
	  	family =  binomial(LNK), ...), 
	  		list(FRM = frm, LNK = lnk))  else
	  substitute(
	  glmer(FRM, data = p.df, 
	  	family =  binomial(LNK), ...),
	  		list(FRM = frm, LNK = lnk))	
	res <- eval(CLL)
#	if (inherits(res, "mer")) print(res, cor = cor) else
#		print(res)
	res
}


