.First.lib <- function(libname, pkgname)
{
	library.dynam("realized")
   
    cat("Realized Library:  Realized Variance, Covariance, Correlation Estimation and Tools.\n")
   # cat("Available in R and S+.\n")
    cat(" R:  http://cran.r-project.org\n")
   # cat("S+:  http://csan.insightful.com\n")
    cat("Copyright (C) 2010  Scott W. Payseur <scott.payseur@gmail.com>\n\n")
   
   # cat("This program is free software with restricted commercial use.\n")
   # cat("See LICENSE file for more details.\n")

   # cat("This program is distributed in the hope that it will be useful,\n")
  #  cat("but WITHOUT ANY WARRANTY; without even the implied warranty of\n")
  #  cat("MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n\n")

    cat("Please report any bugs or feature requests to the author.  User's manual\n")
    cat("in the 'realized' library directory.\n\n\n")
    
    
	if(is.null(version$language)) #SPlus
    {
        data <- function(...)
        {
        }
    }


}


#########################################################################
#
# Utility Functions
#
#########################################################################
.alignedAccum <- function(x,y, period, cum=TRUE, makeReturns...)
{


	x<-.accum.naive(x,x, period)
	y<-.accum.naive(y,y, period)
	
	if(cum)
	{
		ans <- cumsum(x*y)
	}
	else
	{
	    ans <- x*y	
	}
	ans
}


.accum.naive <- function(x,y, period, ...)
{
	.C("rv", 
			as.double(x), #a
			as.double(y), #b
			as.integer(length(x)), #na
			as.integer(period), #period 
			tmpa = as.double(rep(0,as.integer(length(x)/period +1))), #tmp
			as.double(rep(0,as.integer(length(x)/period +1))), #tmp
			as.integer(length(x)/period), #tmpn
			ans = double(1), 
			COPY=c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE), 
			PACKAGE="realized")$tmpa	
}


.alignReturns <- function(x, period, ...)
{
	.C("rv", 
			as.double(x), #a
			as.double(x), #b
			as.integer(length(x)), #na
			as.integer(period), #period 
			tmpa = as.double(rep(0,as.integer(length(x)/period +1))), #tmp
			as.double(rep(0,as.integer(length(x)/period +1))), #tmp
			as.integer(length(x)/period), #tmpn
			ans = double(1), 
			COPY=c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE), 
			PACKAGE="realized")$tmpa	
}

.alignIndices <- function(x, period, ...)
{
	.C("rvperiod", 
			as.double(x), #a
			as.double(x), #b
			as.integer(length(x)), #na
			as.integer(period), #period 
			tmpa = as.double(rep(max(x),as.integer(length(x)/period +1))), #tmp
			as.double(rep(0,as.integer(length(x)/period +1))), #tmp
			as.integer(length(x)/period), #tmpn
			ans = double(1), 
			COPY=c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE), 
			PACKAGE="realized")$tmpa	
}
	
	


#########################################################################
#
# Kernel Estimators 
# See BNHLS (2006), Zhou (1996), HL (1994), HL (1996)
#
#########################################################################


rKernel <- function(x,type=0)
{
	type <- .kernel.chartoint(type)
    .C("justKernel", x=as.double(x),type= as.integer(type), ans=as.double(0),PACKAGE="realized")$ans
}

.kernel.chartoint <- function(type)
{
   if(is.character(type))
    {
    	ans <- switch(casefold(type), 
    	       rectangular=0,
    	       bartlett=1,
    	       second=2,
    	       epanechnikov=3,
    	       cubic=4,
    	       fifth=5,
    	       sixth=6,
    	       seventh=7,
    	       eighth=8,
    	       parzen=9,
    	       th=10,
    	       mth=11,
    	       tukeyhanning=10,
    	       modifiedtukeyhanning=11,
    	       -99)
    	 if(ans==-99)
    	 { 
    	 	warning("Invalid Kernel, using Bartlet")
    	 	1
    	 }
    	 else
    	 {
    	 	ans	
    	 }
    }
    else
    {
    	type
    }
}

rKernel.available <- function()
{
	c("Rectangular", 
	  "Bartlett",
	  "Second",
	  "Epanechnikov",
	  "Cubic",
	  "Fifth",
	  "Sixth",
	  "Seventh",
	  "Eighth",
	  "Parzen",
	  "TukeyHanning",
	  "ModifiedTukeyHanning")
}


rv.kernel <- function(x, q, align.period=1, adj=TRUE, type=0, cts=TRUE, makeReturns=FALSE, rvargs=list(),...)
{	
    if(!is.null(rvargs$period))
    {
    	align.period=rvargs$period
    }
	cdata <- .convertData(x, cts=cts, makeReturns=makeReturns)
	x <- cdata$data
	x <- .alignReturns(x, align.period)
	type <- .kernel.chartoint(type)
	.C("kernelEstimator", as.double(x), as.double(x), as.integer(length(x)),
			 as.integer(q), as.integer(ifelse(adj, 1, 0)),
			 as.integer(type), ab=double(q + 1),
			 ab2=double(q + 1),
			 ans=double(1),PACKAGE="realized")$ans
}

rc.kernel <- function(x, y, q, align.period=1, adj=TRUE, type=0, cts=TRUE, makeReturns=FALSE, ...)
{	
	cdata <- .convertData(x, cts=cts, makeReturns=makeReturns)
	x <- cdata$data
	x <- .alignReturns(x, align.period)
	cdatay <- .convertData(y, cts=cts, makeReturns=makeReturns)
	y <- cdatay$data
	y <- .alignReturns(y, align.period)
	type <- .kernel.chartoint(type)
	.C("kernelEstimator", as.double(x), as.double(y), as.integer(length(x)),
			 as.integer(q), as.integer(ifelse(adj, 1, 0)),
			 as.integer(type), ab=double(q + 1),
			 ab2=double(q + 1),
			 ans=double(1),PACKAGE="realized")$ans
}


#########################################################################
#
# Subsample based estimators
# See AMZ (),(),(), 
#
#########################################################################

.rv.subsample <- function(x, period, cts=TRUE, makeReturns=FALSE,...)
{
	cdata <- .convertData(x, cts=cts, makeReturns=makeReturns)
	x <- cdata$data

	.C("subsample", 
	
			as.double(x), #a
			as.double(x), #na
            as.integer(length(x)), #na
			as.integer(length(x)/period),  	#m
			as.integer(period), #period 
			as.double(rep(0,as.integer(length(x)/period +1))), #tmp
			as.double(rep(0,as.integer(length(x)/period +1))), #tmp
			as.integer(length(x)/period), #tmpn
			ans = double(period), 
			COPY=c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE), 
			PACKAGE="realized")$ans
}


.rc.subsample <- function(x, y, period, cts=TRUE, makeReturns=FALSE, ... )
{
	cdata <- .convertData(x, cts=cts, makeReturns=makeReturns)
	x <- cdata$data

	cdatay <- .convertData(y, cts=cts, makeReturns=makeReturns)
	y <- cdatay$data

	.C("subsample", 
			as.double(x), #a
			as.double(y), #na
            as.integer(length(x)), #na
			as.integer(length(x)/period),  	#m
			as.integer(period), #period 
			as.double(rep(0,as.integer(length(x)/period +1))), #tmp
			as.double(rep(0,as.integer(length(x)/period +1))), #tmp
			as.integer(length(x)/period), #tmpn
			ans = double(period), 
			COPY=c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE), 
			PACKAGE="realized")$ans			
}

rv.timescale <- function(x, period, align.period=1,adj.type="classic", cts=TRUE, makeReturns=FALSE, ...)
{

	x<- .alignReturns(.convertData(x, cts=cts, makeReturns=makeReturns)$data, align.period)

 	n <- dim(as.matrix(x))[[1]]
	nbar <- (n-period+1)/(period)
 	adj <- switch(adj.type, classic=1, adj=(1-(nbar/n))^-1, aa= n/((period-1)*nbar))
    adj * (mean(rv.avg(x, period)) - ((nbar/n) * rv.naive(x,1)))
}

rc.timescale <- function(x,y, period,align.period=1, adj.type="classic", cts=TRUE, makeReturns=FALSE,...)
{
	x<- .alignReturns(.convertData(x, cts=cts, makeReturns=makeReturns)$data, align.period)
	y<- .alignReturns(.convertData(y, cts=cts, makeReturns=makeReturns)$data, align.period)

 	n <- dim(as.matrix(x))[[1]]
	nbar <- (n-period+1)/(period)
 	adj <- switch(adj.type, classic=1, adj=(1-(nbar/n))^-1, aa= n/((period-1)*nbar))
    adj * (mean(rc.avg(x,y, period)) - ((nbar/n) * rc.naive(x,y,1)))
}

rv.avg <- function(x, period, align.period=1, cts=TRUE, makeReturns=FALSE, ...)
{
	x<- .alignReturns(.convertData(x, cts=cts, makeReturns=makeReturns)$data, align.period)
	mean(.rv.subsample(x, period, ...))
}

rc.avg <- function(x, y,  period, align.period=1, cts=TRUE, makeReturns=FALSE, ...)
{
	x<- .alignReturns(.convertData(x, cts=cts, makeReturns=makeReturns)$data, align.period)
	y<- .alignReturns(.convertData(y, cts=cts, makeReturns=makeReturns)$data, align.period)
	mean(.rc.subsample(x, y, period))
}

#########################################################################
#
# Naive estimators
# See ABDL, BNS, etc 
#
#########################################################################
rv.naive <- function(x, period, align.period=1, cts=TRUE, makeReturns=FALSE, ...)
{
	x<- .alignReturns(.convertData(x, cts=cts, makeReturns=makeReturns)$data, align.period)
	.C("rv", 
			as.double(x), #a
			as.double(x), #b
			as.integer(length(x)), #na
			as.integer(period), #period 
			as.double(rep(0,as.integer(length(x)/period +1))), #tmp
			as.double(rep(0,as.integer(length(x)/period +1))), #tmp
			as.integer(length(x)/period), #tmpn
			ans = double(1), 
			COPY=c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE), 
			PACKAGE="realized")$ans	
}

rc.naive <- function(x, y,  period, align.period=1, cts=TRUE, makeReturns=FALSE, ...)
{
	x<- .alignReturns(.convertData(x, cts=cts, makeReturns=makeReturns)$data, align.period)
	y<- .alignReturns(.convertData(y, cts=cts, makeReturns=makeReturns)$data, align.period)

	.C("rv", 
			as.double(x), #a
			as.double(y), #b
			as.integer(length(x)), #na
			as.integer(period), #period 
			as.double(rep(0,as.integer(length(x)/period +1))), #tmp
			as.double(rep(0,as.integer(length(x)/period +1))), #tmp
			as.integer(length(x)/period), #tmpn
			ans = double(1), 
			COPY=c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE), 
			PACKAGE="realized")$ans	
}


#########################################
#
# Global Call 
#
#######################################
rRealized.variance <- function(x, y=NULL, type="naive", period = 1, lags = 1, cor = FALSE, rvargs = list(), cts=TRUE, makeReturns=FALSE)
{
	warning("Depricated:  call rRealizedVarianc")
	rRealizedVariance(x=x, y=y, type=type, period=period, lags=lags, cor=cor, rvargs=rvargs, cts=cts, makeReturns=makeReturns)
}

rRealisedVariance <- function(x, y=NULL, type="naive", period = 1, lags = 1, cor = FALSE, rvargs = list(), cts=TRUE, makeReturns=FALSE)
{
	rRealizedVariance(x=x, y=y, type=type, period=period, lags=lags, cor=cor, rvargs=rvargs, cts=cts, makeReturns=makeReturns)
}


rRealizedVariance <- function(x, y=NULL, type="naive", period = 1, lags = 1, cor = FALSE, rvargs = list(), cts=TRUE, makeReturns=FALSE)
{
	if(!is.null(rvargs$align.period))
	{
	    align.period =rvargs$align.period
	}
	else
	{
	    align.period = 1
	}
	    
#	x<- .alignReturns(.convertData(x, cts)$data, align.period)
#	y<- .alignReturns(.convertData(y, cts)$data, align.period)
	if((n <- .getDataColNum(x)) > 1)
	{
	
	     ans <- matrix(NA, nrow=n, ncol=n)
	     for(i in 1:n)
	     {
	     	for(j in 1:n)
	     	{
	     		if(i == j)
	     		{
	     			if(cor)
	     			{
	     				ans[i,j] <- 1
	     			}
	     			else
	     			{
	     				ans[i,j] <- .realized.variance(x=.getDataCol(x,j), y=NULL, type=type, period=period, lags=lags, rvargs=rvargs, cts=cts, makeReturns=makeReturns)
	     			}
	     		}
	     		else
	     		{
	     		    if(j > i)
	     		    {
	     		    	ans[i,j] <- .realized.variance(x=.getDataCol(x,i), y=.getDataCol(x,j), type=type, period=period, lags=lags, rvargs=rvargs, cor=cor, cts=cts, makeReturns=makeReturns)
	     		    	ans[j,i]<-ans[i,j]
	     		    }	
	     		}
	     	}
	    }
	    ans
	}
	else
	{
		ans <- .realized.variance(x=x, y=y, type=type, period = period, lags = lags, cor = cor, rvargs = rvargs, cts=cts, makeReturns=makeReturns)
	}
#	class(ans) <- "rRealized.variance"
	ans
}



.realized.variance <- function(x, y=NULL, type="naive", period = 1, lags = 1, cor = FALSE, rvargs = list(), cts=TRUE,makeReturns=FALSE)
{   
	if(cor)
	{
		rvx <- do.call(paste("rv.", type, sep=""), c(rvargs=rvargs,list(x=x, q=lags, k=lags, period=period, cts=cts, makeReturns=makeReturns)))
		rvy <- do.call(paste("rv.", type, sep=""), c(rvargs=rvargs,list(x=y, q=lags, k=lags, period=period, cts=cts, makeReturns=makeReturns)))
		rcxy <- do.call(paste("rc.", type, sep=""),c(rvargs=rvargs,list(x=x, y=y, q=lags, k=lags, period=period, cts=cts, makeReturns=makeReturns)))
		rcxy/(sqrt(rvx)*sqrt(rvy))
	}
	else
	{
		funct <- paste(ifelse(is.null(y), "rv.", "rc."), type, sep="")
		do.call(funct, c(rvargs=rvargs, list(x=x, y=y, q=lags, k=lags, period=period, cts=cts, makeReturns=makeReturns)))	
	}
}

rc.zero <- function(x, y, period, align.period=1, cts=TRUE, makeReturns=FALSE, ...)
{
	y<- .alignReturns(.convertData(y, cts=cts, makeReturns=makeReturns)$data, align.period)
	x<- .alignReturns(.convertData(x, cts=cts, makeReturns=makeReturns)$data, align.period)
	acy <- .accum.naive(x=y,y=y,period=period)
	acx <- .accum.naive(x=x,y=x,period=period)
    sum((acx*acy)==0)/length(acy)
}

rv.zero <- function(x, period, align.period=1, cts=TRUE, makeReturns=FALSE, ...)
{
	x<- .alignReturns(.convertData(x, cts=cts, makeReturns=makeReturns)$data, align.period)
	ac <- .accum.naive(x=x,y=x,period=period)
    sum((ac*ac)==0)/length(ac)
}



rCumSum <- function(x, period = 1, align.period=1, plotit=FALSE, type='l', cts = TRUE, makeReturns=FALSE)
{
	
	ans <- list(x = NULL, y = NULL)
    ans$x <- .alignIndices(1:length(.convertData(x, cts=cts, makeReturns=makeReturns)$data), align.period)
    ans$x <- .alignIndices(ans$x, period)

	x<- .alignReturns(.convertData(x, cts=cts, makeReturns=makeReturns)$data, align.period)
	x<- .alignReturns(.convertData(x, cts=cts, makeReturns=makeReturns)$data, period)
	
	ans$y <- cumsum(x)
	if(plotit)
	{
		plot(cumsum(x), xlab="Time", ylab="Cummulative Returns", type=type)
		return(NULL)
	}
	ans
}

rScatterReturns <- function(x,y, period, align.period=1,numbers=FALSE,xlim= NULL, ylim=NULL, plotit=TRUE, pch=NULL, cts=TRUE, makeReturns=FALSE, scale.size=0, col.change=FALSE,...)
{

	y<- .alignReturns(.convertData(y, cts=cts, makeReturns=makeReturns)$data, align.period)
	x<- .alignReturns(.convertData(x, cts=cts, makeReturns=makeReturns)$data, align.period)
	
    x<-.accum.naive(x, x, period)
	y<-.accum.naive(y, y, period)
	if(is.null(pch))
	    pch=1
	
	it <- table(round(x,4),round(y,4))
	xs <- as.numeric(dimnames(it)[[1]])
	ys <- as.numeric(dimnames(it)[[2]])

	if(is.null(ylim))
	    ylim=c(min(ys), max(ys))
	if(is.null(xlim))
	    xlim=c(min(xs), max(xs))
	
	mat <- matrix(it, nrow=length(xs), ncol=length(ys))
	
	if(plotit)
	{
		plot(0,0, xlim=xlim, ylim=ylim , type='n',...)
		lines(c(0,0), c(-1,2), col="grey", lty=3, lwd=2)
		lines(c(-1,2), c(0,0), col="grey", lty=3, lwd=2)
	
        	maxed <- max(mat)

		for(i in 1:length(xs))
		{
    		for(j in 1:length(ys))
    		{
            	if(mat[i,j]!=0)
            	{
            		if(col.change)
            		   thecol <- round(runif(1)*100,0)
            		else
            		   thecol = 1
            		   
        			if(numbers)
        			{
				     
				    	if(scale.size ==0)
							text(xs[i], ys[j],as.character(mat[i,j]), cex=.7, col=thecol)    	
                        else
                        	text(xs[i], ys[j], as.character(mat[i,j]), cex = (mat[i,j]/maxed) * scale.size, col=thecol)
        			}
        			else
        			{
        				if(scale.size ==0)
                			points(xs[i], ys[j], pch=pch, cex=.7, col=thecol)    	
                        else
                			points(xs[i], ys[j], pch=pch, cex = (mat[i,j]/maxed) * scale.size, col=thecol)
        			}
            	}

    		}
		}
		return(NULL)
	
	}	
	mat
}
    




rc.hy <- function(x,y, period=1, align.period =1, cts = TRUE, makeReturns=FALSE, ...)
{
	cdata <- .convertData(x, cts=cts, makeReturns=makeReturns)
	x <- cdata$data
	x.t <- cdata$milliseconds

	cdatay <- .convertData(y, cts=cts, makeReturns=makeReturns)
	y <- cdatay$data
	y.t <- cdatay$milliseconds
	
	
	errorCheck <- c(is.null(x.t),is.na(x.t), is.null(y.t), is.na(y.t))
	if(any(errorCheck))
	    stop("ERROR: Time data is not in x or y.")
	    

    sum(	.C("pcovcc", 
			as.double(x), #a
			as.double(rep(0,length(x)/(period*align.period)+1)),
			as.double(y), #b
			as.double(x.t), #a
			as.double(rep(0,length(x)/(period*align.period)+1)), #a
			as.double(y.t), #b
			as.integer(length(x)), #na
			as.integer(length(x)/(period*align.period)),
			as.integer(length(y)), #na
			as.integer(period*align.period),
			ans = double(length(x)/(period*align.period)+1), 
			COPY=c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE), 
			PACKAGE="realized")$ans)
}


#######################################################################
#
# Graphs 
#
#######################################################################

#
# Signature Plot
#
rSignature <- function(range, x, y=NULL, type="naive", cor = FALSE, rvargs = list(), xscale=1, iteration.funct="", iterations=NULL, plotit=FALSE, cts=TRUE, makeReturns=FALSE)
{
	
    
	if(is.null(iterations))
	{
		x <- .convertData(x, cts=cts, makeReturns=makeReturns)$data
        y <- .convertData(y, cts=cts, makeReturns=makeReturns)$data
    	ans <- list(x=range*xscale, 
        	 y=sapply(range, function(r, x, y, type, period, lags, cor, rvargs){.realized.variance(x=x, y=y, type=type, period = r, lags = r, cor = cor, rvargs = rvargs)},x=x,y=y,type=type, cor=cor, rvargs=rvargs),
        	 xgrid=range,
        	 type = type,
        	 cor = cor,
        	 cov = is.null(y),
        	 cts= cts)	
	}
	else
	{
	     ans <- sapply(1:length(iterations), 
	     	function(i,iterations, iteration.funct,range, x, y, type, period, lags, cor, rvargs, cts, makeReturns)
	     	{
	     		if(!is.null(y))
	     		{
	     	    	y.new <- do.call(iteration.funct, args=list(x=y, i=iterations[i], rvargs=rvargs))
	     		}
	     		else
	     		{
	     	    	y.new <- NULL
	        	}
	     		sapply(range, 
	     	 		function(r, x, y, type, period, lags, cor, rvargs, cts, makeReturns)
	         		{
	     	     		.realized.variance(x=x, y=y, type=type, period = r, lags = r, cor = cor, rvargs = rvargs)
	         		},
	         	x=do.call(iteration.funct, args=list(x=x, i=iterations[i], rvargs=rvargs)), 
	         	y=y.new,
	         	type=type, 
	         	cor=cor, 
	         	rvargs=rvargs, 
	         	cts=cts,
	         	makeReturns=makeReturns)
	     	},
	     	 iterations=iterations, 
	     	 range=range, 
	     	 iteration.funct=iteration.funct,
	     	 x=x,
	     	 y=y,
	     	 type=type, 
	     	 cor=cor, 
	     	 rvargs=rvargs, 
	     	 cts=cts,
	     	 makeReturns=makeReturns)
             if(class(ans)!="numeric")
             {
             	ans <- apply(ans,1,mean)
             }
    	ans <- list(x=range*xscale, 
        	 y=ans,
        	 xgrid=range,
        	 type = type,
        	 cor = cor,
        	 cov = is.null(y))	

    }
    #class(ans) <- "rSignature"
    if(plotit)
    {
        .rSignature.plot(ans)
    }
    ans
}

.rSignature.plot <- function(obj)
{
    if(obj$cov && obj$cor)
    {
        ylab = "Realized Correlation"
    }
    else{
    	if(obj$cov)
    	    ylab = "Realized Covariance"
    	else
    	    ylab = "Realized Variance"
    }
    xlab = "Sampling Frequency"
    main = paste(ylab, ":", obj$type, sep="")
    plot(obj$x, obj$y, xlab=xlab, ylab=ylab, main=main)		
}


#
# Helper functions
#
.getDataColNum <- function(x)
{
	
	if(class(x) == "realizedObject")
    {
        return(dim(x)[[2]])	
    }
    if(is.null(version$language)) #splus
    {
	     if(class(x) == "timeSeries")
	     {
	         return(dim(x)[[2]])
	     }
    }

    if(class(x) == "list")
    {
    	if(is.null(x$data))
    	{
    		return(NA)
    	}
    	else
    	{
    		if(class(x$data)=="matrix" || class(x$data)=="data.frame")
    		{
    			return(dim(x$data)[[2]])
    		}
    		else
    		{
    			return(1)
    		}
    	}
    }
    else
    {	
    	if(class(x)=="matrix" || class(x)=="data.frame")
    	{
    		return(dim(x)[[2]])
    	}
    	else
    	{
    		return(1)
    	}
    }
}


.getDataCol <- function(x,i)
{
	if(is.null(x))
	{
		return(NULL)
	}


    if(class(x) == "realizedObject")
    {
        return(x[,i])	
    }	
	if(is.null(version$language)) #splus
    {
	     if(class(x) == "timeSeries")
	     {
	         return(x[,i])
	     }
    }
	
    if(class(x) == "list")
    {
    	if(is.null(x$data))
    	{
    		return(x)
    	}
    	else
    	{
    		if(class(x$data)=="matrix" || class(x$data)=="data.frame")
    		{
    			return(x$data[,i])
    		}
    		else
    		{
    			return(x$data)
    		}
    	}
    }
    else
    {	
    	if(class(x)=="matrix" || class(x)=="data.frame")
    	{
    		return(x[,i])
    	}
    	else
    	{
    		return(x)
    	}
    }
}

rMarginal <- function(x, y=NULL, period, align.period=1, plotit=FALSE, cts=TRUE, makeReturns=TRUE)
{
	ans <- list(x = NULL, y = NULL)
    ans$x <- .alignIndices(1:length(x), align.period)
    ans$x <- .alignIndices(ans$x, period)

	if(is.null(y))
	    y <- x
	    
	x<- .alignReturns(.convertData(x, cts=cts, makeReturns=makeReturns)$data, align.period)
	y<- .alignReturns(.convertData(y, cts=cts, makeReturns=makeReturns)$data, align.period)
	

    ans$y <- .alignedAccum(x=x, y=y, period=period, cum=FALSE)
    
    if(plotit)
    {
    	plot(ans, xlab="", ylab="Realized Marginals")
        return(NULL)
    }
    ans
}

rAccumulation <- function(x, period=1, y=NULL, align.period=1, plotit=FALSE, cts=TRUE, makeReturns=FALSE)
{
    ans <- list(x=NULL, y=NULL)
    ans$y <- cumsum(rMarginal(x=x, y=y, period=period, align.period=align.period, cts=cts, makeReturns=makeReturns)$y)
#    ans$x <- .alignIndices(1:length(x), align.period)
 #   ans$x <- .alignIndices(ans$x, period)
    ans$x <- rCumSum(x=x, period=period, align.period=align.period, cts=cts, makeReturns=makeReturns)$x
    #ans$x <- ans$x[-length(ans$x)]
    if(plotit)
    {
    	plot(ans, xlab="", ylab="Realized Accumulation")
    	return(NULL)
    }
    ans
}
 



realizedObject <- function(x, cts = TRUE, makeReturns = FALSE, millisstart=34200000, millisend=57600000)
{
	
    y <- .convertData(x, cts=cts, millisstart=millisstart, millisend=millisend, makeReturns=makeReturns)	
    class(y) <- "realizedObject"
    y$cts = cts
#    if(makeReturns)
#    {
#    	errcheck <- try(getReturns(y$data))
#        if(class(errcheck) != "Error")
#        {
#        	y$data <- errcheck
#        }
#        else
#        {
#        	warning("It appears that these are already returns.  Not creating returns")
#        }
#    }
    y
}

print.realizedObject <- function(x, ...)
{
    if(class(x$data) == "matrix")
    {
    	n <- dim(x$data)[[1]]
    	k <- dim(x$data)[[2]]
        display.n <- ifelse(n < 10, n, 10)
        display.obj <- cbind(x$data[1:display.n,], x$milliseconds[1:display.n])
        dims <- dimnames(x$data)[[2]]
        if(is.null(dims))
            dims <- paste("data", 1:k, sep="")
           
	    dimnames(display.obj) <- list(rep("",display.n), c(dims, "milliseconds"))
    }
    else
    {
        n <- length(x$data)
        k <- 1
        display.n <- ifelse(n < 10, n, 10)
        display.obj <- matrix(c(x$data[1:display.n], x$milliseconds[1:display.n]), ncol=1+k)
    	dimnames(display.obj) <- list(rep("",display.n), c("data", "milliseconds"))
	
    }
    cat("Realized Object: (length = ", n, ", cts=", x$cts,")\n")
	print(display.obj)
	cat("...")
	
	
}


"[.realizedObject" <- function(x, i=NULL,j=NULL, drop=TRUE)
{
	ret <- x
	ret$milliseconds <- x$milliseconds[i]
	if(class(x$data) == "matrix")
    {
        if(is.null(i))
            i <- 1:(dim(x$data)[[1]]) 
        if(is.null(j))
            j <- 1:(dim(x$data)[[2]]) 
		ret$data <- x$data[i,j]
    }
    else
    	ret$data <- x$data[i]
	ret
}

dim.realizedObject <- function(x)
{
    if(class(x$data) == "matrix")
        return(dim(x$data))
    else
        return(c(length(x$data), 1))
}
merge.realizedObject <- function(x, y=NULL,...)
{
	if(is.null(y))
	{
		return(x)
	}

    k <- nargs()
    inputs <- list(x,y, ...)
    inputs.class <- sapply(1:k, function(x, inputs){class(inputs[[x]])}, inputs)
	if(sum(inputs.class != "realizedObject"))
	{
		stop("merge.realizedObject takes object of type realizedObject only.")
	}
	inputs.len <- sapply(1:k, function(x, inputs){length(inputs[[x]]$milliseconds)}, inputs)
	if(sum(inputs.len != inputs.len[[1]]))
	{
		stop("Cannot merge objects with different timings")
	}
	if(k == 1)
	{
		return(x)
	}
	else
	{
		for(i in 2:k)
		{
		    x$data <- cbind(x$data,inputs[[i]]$data)
		}
	}
	x
}

#
# Data Handling
#
.convertData <- function(x, cts = TRUE, millisstart=34200000, millisend=57600000, makeReturns=FALSE)
{
	if(is.null(x))
	{
		return(NULL)
	}
	if(class(x) == "realizedObject")
	{
		return(x)
	}
    if(is.null(version$language)) #splus
    {
	     if(class(x) == "timeSeries")
	     {
	     	x <- x[!is.na(x[,1]),1]
	     	if(cts)
	     	{
	     	    return(ts2realized(x, millisstart=millisstart, millisend=millisend, make.returns=makeReturns)$cts)
	     	}
	     	else
	     	{
	     	    return(ts2realized(x, millisstart=millisstart, millisend=millisend, make.returns=makeReturns)$tts)
	     	}
	     	#list(milliseconds = positions(x)@.Data[[2]], data = matrix(seriesData(x), ncol=1))
	     }
    }
	
	if(class(x) == "list")
	{
	 	if(sum(names(x) == c("tts", "cts")) == 2) #realized obj  
		{
		    if(cts)
	 	    {
	           return(x$cts)
	        }
	     	else
	     	{
	            return(x$tts)
            }
         }
         if(sum(names(x) == c("data", "milliseconds")) == 2) # realized object cts or tts, just return
         {
       	    if(makeReturns)
   		    {
    		  	errcheck <- try(getReturns(.sameTime(x$data, x$milliseconds)))
       			if(class(errcheck) != "Error")
       			{
       				x$data <- errcheck
       		    	x$milliseconds <- intersect(x$milliseconds,x$milliseconds)
       			}
		      	else
       			{
       				warning("It appears that these are already returns.  Not creating returns")
       			}
    	    }		
  		    else
  		    {
    	    	x$data <- .sameTime(x$data, x$milliseconds)
    	    	x$milliseconds <- intersect(x$milliseconds,x$milliseconds)
   		    }		
   	    	if(cts)
        	{
     		    toret <- list(data=.toCts(x=x$data, millis=intersect(x$milliseconds,x$milliseconds), millisstart=millisstart, millisend=millisend),
	                          milliseconds=(((millisstart/1000)+1):(millisend/1000))*1000)
	            return(toret)
	     	}
	     	else
	     	{
	     	    toret <- list(data=x$data, 
	     	                  milliseconds=intersect(x$milliseconds,x$milliseconds))
	     	    return(toret)
	     	 }
	     }
	}
	
	
	if(class(x) == "timeSeries")
    {
        stop("R timeSeries not implmented yet. Convert to realized object")
    }
    return(list(milliseconds = 1:dim(as.matrix(x))[[1]], data = as.matrix(x)))  # not an object, fake the milliseconds and return
}

plot.realizedObject <- function(x,y=NULL,...)
{
    plot(x$milliseconds, x$data, xlab="Milliseconds", ylab="", main="")
    NULL
}


getReturns <- function(x)
{
   	   n <- length(x)[[1]]
   	   return(log(x[2:n]) - log(x[1:(n-1)]))
}




timeDate <- function(x, format)
{
   warning("This function is for SPLUS and does not work")
   x
}



.ts2millis <- function(x,...)
{
	
	millis <- as.numeric(as.character((timeDate(as.character(x@positions), format="%H")))) * 60 * 60 * 1000 +
	          as.numeric(as.character((timeDate(as.character(x@positions), format="%M")))) * 60 * 1000 +
	          as.numeric(as.character((timeDate(as.character(x@positions), format="%S")))) * 1000 +
              as.numeric(as.character((timeDate(as.character(x@positions), format="%N"))))
	millis
}
		
		
ts2realized <- function(x, make.returns=TRUE,millisstart=34200000, millisend=57600000)
{
	thedata <- data.sameTime(as.numeric(as.matrix(x@data)), .ts2millis(x))

    if(make.returns)
    {
    	
		thedata <- getReturns(thedata)
    	
		tts <- list(data=as.numeric(thedata), milliseconds=intersect(.ts2millis(x),.ts2millis(x))[-1])
		cts <- list(data=.toCts(x=as.numeric(thedata), millis=intersect(.ts2millis(x),.ts2millis(x)), millisstart=millisstart, millisend=millisend),
	    	        milliseconds=(((millisstart/1000)+1):(millisend/1000))*1000)
    }
    else
    {
		tts <- list(data=as.numeric(thedata), milliseconds=intersect(.ts2millis(x),.ts2millis(x)))
		cts <- list(data=.toCts(x=as.numeric(thedata), millis=intersect(.ts2millis(x),.ts2millis(x)), millisstart=millisstart, millisend=millisend),
	    	        milliseconds=(((millisstart/1000)+1):(millisend/1000))*1000)
    	
    	
    }
    	ans <- list(tts=tts, cts=cts)	
	ans
}

		
tsGetDay<-function(ts, dateString)
{
	if(is(ts, "timeSeries"))
		pos = ts@positions
	else stop("ts must be a timeSeries object")
	pos@format = "%02m/%02d/%Y"
	poschar = as.character(pos)
	inds <- poschar == dateString
	ts[inds]
}
		
		
tsGetDayObject<-function(x, i, cts=TRUE, ...)
{
	ts = x
	dateString=i
	if(cts)
    ts2realized(tsGetDay(ts, dateString))$cts$data
    else
    ts2realized(tsGetDay(ts, dateString))$tts$data    
}


#tsGetDayObject(msftt.ts[,6], "05/01/1997")
		#modeule(finmetrics)
#		yoyo<-(tsGetDay(msftt.ts, "05/01/1997")[,6])
   #     theday <- tsGetDay(msftt.ts, "05/01/1997")[,6]
    #    yoyo <- data.sameTime(as.numeric(as.matrix(theday@data)), .ts2millis(theday))
     #   yoyo <- getReturns(yoyo)
        
        
		
data.sameTime <- function(x, millis)
{
	.sameTime(x=x,millis=millis)
}

.sameTime <- function(x, millis)
{
       .C("sametime", 
			as.double(x), #a
			as.integer(length(x)), #na
			as.integer(millis), #millis
			ans = double(length(union(millis,millis))), #tts
			COPY=c(FALSE,FALSE,FALSE,TRUE), 
			PACKAGE="realized")$ans
}




data.toCts <- function(x, millis, millisstart=34200000, millisend=57600000)
{
	.toCts(x=x, millis=millis, millisstart=millisstart, millisend=millisend)
}

.toCts <- function(x, millis, millisstart=34200000, millisend=57600000)
{
       .C("tocts", 
			as.double(x), #a
			as.integer(length(x)),
			as.integer(millis), #millis
            as.integer(millisstart),
            as.integer(millisend),
			ans = double(((millisend-millisstart)/1000)), #cts
			COPY=c(FALSE,FALSE,FALSE,FALSE,TRUE), 
			PACKAGE="realized")$ans
}

data.toReturns <- function(x)
{
    x <- as.numeric(x)   
    n <- length(x)
    log(x[2:n]) - log(x[1:(n-1)])
}











#########################################################################
#
# Optimal Timing
# 
#
#########################################################################

#
#rv.star <- function(x,period,...)
#{
#    rv.naive(x,period=period.star(x,quarticity.period=period))	
#}
#
#rc.star <- function(x, y, period,...)
#{
#	print(period.star(x,y=y,quarticity.period=period))
#    rc.naive(x,y=y,period=period.star(x,y=y,quarticity.period=period,exact=FALSE))	
#}
#
#noise.fourth <- function(x, period=1, align.period=1, y = NULL, lead = FALSE, ...)
#{
#	cdata <- .convertData(x)
#	x <- cdata$data
#	x <- .alignReturns(x, align.period)
#    if(is.null(y))
#    {
#        y <- x
#    }
#    else
#    {
#        cdatay <- .convertData(y)
#	    y <- cdatay$data
#	    y <- .alignReturns(y, align.period)	
#    }    
#    if(lead)
#    {
#       .C("rfourthlead", 
#			as.double(x), #a
#			as.double(y), #b
#			as.integer(length(x)), #na
#			as.integer(period), #period 
#			as.double(rep(0,as.integer(length(x)/period +1))), #tmp
#			as.double(rep(0,as.integer(length(x)/period +1))), #tmp
#			as.integer(length(x)/period), #tmpn
#			ans = double(1), 
#			COPY=c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE), 
#			PACKAGE="realized")$ans	
#    	
#    }
#    else
#    {
#	    .C("rfourth", 
#			as.double(x), #a
#			as.double(y), #b
#			as.integer(length(x)), #na
#			as.integer(period), #period 
#			as.double(rep(0,as.integer(length(x)/period +1))), #tmp
#			as.double(rep(0,as.integer(length(x)/period +1))), #tmp
#			as.integer(length(x)/period), #tmpn
#			ans = double(1), 
#			COPY=c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE), 
#			PACKAGE="realized")$ans	
#    }
#}
#
#noise.second <- function(x, period=1, align.period=1, y = NULL,...)
#{
#   cdata <- .convertData(x)
#   x <- cdata$data
#   x <- .alignReturns(x, align.period)
#       
#   if(is.null(y))
#   {
#        y <- x
#   }
#   else
#   {
#       cdatay <- .convertData(y)
#       y <- cdatay$data
#       y <- .alignReturns(y, align.period)
#   }
#   rc.naive(x=x,y=y,period=period)/length(x)
#}
#
#
#quarticity <- function(x, period = 900, align.period=1, y=NULL)
#{
#	cdata <- .convertData(x)
#    x <- cdata$data
#    x <- .alignReturns(x, align.period)
#    M <- length(x)/period
#	if(is.null(y))
#	{
#	    (M/3) * noise.fourth(x=x, period=period)    	
#	}
#	else
#	{
#    	cdatay <- .convertData(y)
#        y <- cdatay$data
#        y <- .alignReturns(y, align.period)
#
#	    mbar <- length(x)/period
#	    mbar * noise.fourth(x=x, period=period, y=y, lead=FALSE) -  mbar * noise.fourth(x=x, period=period, y=y, lead=TRUE)
#	}
#}
#
#MSE <- function(x, period, y = NULL, align.period=1, quarticity.period=900, quarticity.x=x, quarticity.y=y, M=NULL)
#{
#    cdata <- .convertData(x)
#    x <- cdata$data
#    x <- .alignReturns(x, align.period)
#
#    cdata.q <- .convertData(quarticity.x)
#    quarticity.x <- cdata.q$data
#    quarticity.x <- .alignReturns(quarticity.x, align.period)
# 
# 
#    M <- length(x)/ period	
#    	
#	if(!is.null(y))
#	{
#        cdatay <- .convertData(y)
#	    y <- cdatay$data	
#        y <- .alignReturns(y, align.period)	    
#
#        cdata.qy <- .convertData(quarticity.y)
#        quarticity.y <- cdata.qy$data
#        quarticity.y <- .alignReturns(quarticity.y, align.period)
#        
#        
#    	M * noise.fourth(x=x, period=1, y=y) + 2*(M-1)* noise.fourth(x=x, period=1, y=y, lead=TRUE) 
#    	+ (M^2 - 3*M + 2) * noise.second(x=x, period=1, y=y)^2 + noise.second(x=x, period=1) * rv.naive(x=x, period=period) +
#    	2 * noise.second(x=x, period=1, y=y) * rc.naive(x=x, period=period, y=y) +noise.second(x=y, period=1) * rv.naive(x=y, period=period) +
#    	quarticity(x=quarticity.x, y=quarticity.y, period=quarticity.period)/M	    
#	}
#	
#	else
#	{
#    	B <- 2 * noise.fourth(x, period=1) - 3 * noise.second(x, period=1)^2
#    	a <- noise.second(x, period=1)^2
#    	g <- (4 * noise.second(x, period=1) * rv.naive(x, period=period)) - noise.fourth(x,period=1) +
#        	  2*noise.second(x, period=1)^2	
#		(2 * quarticity(quarticity.x, period = quarticity.period))/M + M*B +M*M*a + g  
#	}
#}
#
#period.star <- function(x, y = NULL, align.period=1, quarticity.period=900, interval=c(1,1000), exact=TRUE,quarticity.x=x, quarticity.y=y,...)
#{
#
#	    cdata <- .convertData(x)
#    	x <- cdata$data
#        x <- .alignReturns(x, align.period)
#
#        cdata.q <- .convertData(quarticity.x)
#        quarticity.x <- cdata.q$data
#        quarticity.x <- .alignReturns(quarticity.x, align.period)
#
#
#		if(!is.null(y))
#		{
#			cdatay <- .convertData(y)
#        	y <- cdatay$data
#        	y <- .alignReturns(y, align.period)
#		
#		    cdata.qy <- .convertData(quarticity.y)
#            quarticity.y <- cdata.qy$data
#            quarticity.y <- .alignReturns(quarticity.y, align.period)
#       
#		}
#    	if(exact)
#	{
#		if(is.null(y))
#		{
#		
#		     round(optimize(function(p, x, quarticity.period, quarticity.x){MSE(x,round(p,0),quarticity.period=quarticity.period, quarticity.x=quarticity.x)},
#	    	                              interval=interval, 
#	        	                          x=x, 
#	            	                      quarticity.period=quarticity.period,
#	            	                      quarticity.x=quarticity.x
#	            	                  )$minimum)
#		}
#		else
#		{
#		    round(optimize(function(p, x, y, quarticity.period, quarticity.x, quarticity.y){MSE(x,period=round(p,0),y=y,quarticity.period=quarticity.period,quarticity.x=quarticity.x,quarticity.y=quarticity.y)},
#	    	                              interval=interval, 
#	        	                          x=x,
#	        	                          y=y, 
#	            	                      quarticity.period=quarticity.period,
#	            	                      quarticity.x=quarticity.x,
#	            	                      quarticity.y=quarticity.y)$minimum,...)	
#			
#		}
#	}
#	else
#	{
#    	
#
##   		M <- length(x)/ period	
#
#   		
#   		if(is.null(y))
#		{
#			print(noise.fourth(quarticity.x, period=quarticity.period))
#			print(noise.second(x,period=1))
#   			round(length(x)/(((1/3)*noise.fourth(quarticity.x, period=quarticity.period)) / noise.second(x,period=1)^2)^(1/3))
#		}
#		else
#		{
#			round(length(x)/ (quarticity(x=quarticity.x, y=quarticity.y, period=quarticity.period)/(2*noise.second(x=x,y=y,period=1)^2))^(1/3))
#		}
#	}
#}
#
#m.star <- function(x, y = NULL, quarticity.period=900, align.period=1, interval=c(1,1000), exact=TRUE,quarticity.x=x, quarticity.y=y,...)
#{
#    cdata <- .convertData(x)
#    x <- cdata$data
#    x <- .alignReturns(x, align.period)
#    if(!is.null(y))
#    {
#        cdatay <- .convertData(y)
#        y <- cdatay$data
#        y <- .alignReturns(y, align.period)
#    	
#    }  
#    round(length(x)/ period.star(x, y=y, quarticity.period=quarticity.period, interval=interval, exact=exact, quarticity.x=x, quarticity.y=y,...),0)
#    	
#}
#
#phi.star <- function(x,  quarticity.period=900, interval=c(1,1000), exact=TRUE,quarticity.x=x)
#{
#    
#    cdata <- .convertData(x)
#    x <- cdata$data  
#    M <- length(x)  
#    period <- period.star(x,  quarticity.period=quarticity.period, interval=interval, exact=exact)
#    v <- rv.naive(x, period=period)
#    Q <- quarticity(quarticity.x, period=quarticity.period)
#    
#    ( 1.5 * (( v^2/ M^2 ) / Q) ) ^ (1/3)
#}
#
#lag.star <- function(x,  quarticity.period=900, interval=c(1,1000), exact=TRUE, quarticity.x=x, cts=TRUE)
#{
#    
#    cdata <- .convertData(x, cts=cts)
#    x <- cdata$data  
#
# 	floor(phi.star(x, quarticity.period=quarticity.period, interval=interval, exact=exact,quarticity.x=quarticity.x) * length(x))
#}
#
