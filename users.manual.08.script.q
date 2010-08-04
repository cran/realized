#
# # testRealized.R			
#
# author: Scott Payseur
# created: Sept 12, 2007
#
# Data used in examples:
# Data are available for download from
#
#

library(realized)
data(msft.real.cts)
data(msft.real.tts)
data(ge.real.cts)
data(ge.real.tts)

msft.real.cts[[1]][1:10]
msft.real.cts[[1]][1:5] 
 
 
# Figure 1
#java.graph("d:/dev/realizedDoc/Beta08/UsersManual/plots/Figure1.jpg", format="JPEG", width=1400, height=1200)
par(mfrow=c(2,3))
tmp <- sapply(1:6, function(x, rets){plot(rCumSum(rets[[x]]), ylab="Cumulative Return", xlab="")}, rets=msft.real.tts)
#dev.off() 
 
data(dates.example)
dates.example
 
#
# Figure 2
#
par(mfrow=c(1,1))
simpleIteration <- function(x, i,args){x[[i]]}
#java.graph("d:/dev/realizedDoc/Beta08/UsersManual/plots/Figure2.jpg", format="JPEG", width=1400, height=1200)
plot(rSignature((1:120)*10+1, msft.real.cts, xscale=1/60, iteration.funct="simpleIteration", iterations=1:6),
                ylab="Realized Variance", xlab="Sampling Frequency (Minutes)", main="MSFT", sub=paste(dates.example[[1]], dates.example[[6]], sep=" - "))
#dev.off()                 
args(rSignature)

test.sig <- rSignature(1:1200, msft.real.cts[[1]], xscale=1/60)
names(test.sig)

#java.graph("d:/dev/realizedDoc/Beta08/UsersManual/plots/Figure3A.jpg", format="JPEG", width=1400, height=1200)
plot(test.sig, ylab="Realized Variance", xlab="Sampling Frequency (Minutes)", main="MSFT", sub=dates.example[[1]])
#dev.off()
#java.graph("d:/dev/realizedDoc/Beta08/UsersManual/plots/Figure3B.jpg", format="JPEG", width=1400, height=1200)
plot(x=test.sig$x[-(1:20)], y=test.sig$y[-(1:20)], ylab="Realized Variance", xlab="Sampling Frequency (Minutes)", main="MSFT", sub=dates.example[[1]])
#dev.off()

# Figure 4A
#java.graph("d:/dev/realizedDoc/Beta08/UsersManual/plots/Figure4A.jpg", format="JPEG", width=1400, height=1200)
acf(msft.real.cts[[1]]$data, main="ACF: MSFT")
#dev.off()
 
# Figure 4B
#java.graph("d:/dev/realizedDoc/Beta08/UsersManual/plots/Figure4B.jpg", format="JPEG", width=1400, height=1200)
plot(x=test.sig$x[-(1:20)],y=test.sig$y[-(1:20)],ylab="Realized Variance", xlab="Sampling Frequency (Minutes)", main="MSFT",sub=dates.example[[1]])
test.rect <- rSignature(1:400, msft.real.cts[[1]], xscale=1/20, type="kernel", args=list(type="rectangular"))
lines(test.rect, col=2, lwd=2) 
axis(3, c(0,(1:5)*4), c("Lags:",as.character((1:5)*80))) 
legend(15,.0008,c("Rectangular"), lwd=c(2), col=c(2)) 
#dev.off()
rKernel.available()
 
# Figure 5
#java.graph("d:/dev/realizedDoc/Beta08/UsersManual/plots/Figure5.jpg", format="JPEG", width=1400, height=1200)
par(mfrow=c(3,4)) 
x <- (0:100)*.01 
for(i in 1:length(rKernel.available())) 
   plot(x=x,y=sapply(x, FUN="rKernel", type=rKernel.available()[i]), xlab="", ylab="", main=rKernel.available()[i],ylim=c(0,1))
#dev.off()

# Figure 6a
#java.graph("d:/dev/realizedDoc/Beta08/UsersManual/plots/Figure6A.jpg", format="JPEG", width=1400, height=1200)
par(mfrow=c(1,1))
plot(x=test.sig$x[-(1:20)],y=test.sig$y[-(1:20)],ylab="Realized Variance", xlab="Minutes", main="MSFT",sub=dates.example[[1]])
test.mth <- rSignature(1:400, msft.real.cts[[1]], xscale=1/20, type="kernel", args=list(type="mth"))
test.bart <- rSignature(1:400, msft.real.cts[[1]], xscale=1/20, type="kernel", args=list(type="bartlett"))
lines(test.mth, col=3, lwd=2) 
lines(test.bart, col=4, lwd=2)
axis(3, c(0,(1:5)*4), c("Lags:",as.character((1:5)*80))) 
legend(15,.0008,c("Mod T-H", "Bartlett"), lwd=c(2,2), col=c(3,4)) 
#dev.off()

#java.graph("d:/dev/realizedDoc/Beta08/UsersManual/plots/Figure6B.jpg", format="JPEG", width=1400, height=1200)
# Figure 6b
test.sig.min <- rSignature(1:120, msft.real.cts[[1]], xscale=1/2, args=list(align.period=30)) 
plot(test.sig.min, ylab="Realized Variance", xlab="Minutes", main="MSFT", sub=dates.example[[1]])
test.tt <- rSignature(1:20, msft.real.cts[[1]], xscale=3, type="timescale", args=list(adj.type="classic", align.period=60)) 
test.tt.adj <- rSignature(1:20, msft.real.cts[[1]], xscale=3, type="timescale", args=list(adj.type="adj", align.period=60)) 
test.tt.aa <- rSignature(1:20, msft.real.cts[[1]], xscale=3, type="timescale", args=list(adj.type="aa", align.period=60)) 
lines(test.tt, col=3, lwd=2) 
lines(test.tt.adj, col=4, lwd=2) 
lines(test.tt.aa, col=5, lwd=2) 
axis(3, c(0,(1:5)*12), c("Subgrids:",as.character((1:5)*4))) 
legend(45,.0006,c("Classic", "Adj", "AA"), lwd=c(2,2,2), col=c(3,4,5)) 
#dev.off()
 
# Traditional Estimate at highest frequency 
rRealizedVariance(x=msft.real.cts[[1]], type="naive", period=1) 

# Traditional Estimate at one minute frequency 
rRealizedVariance(x=msft.real.cts[[1]], type="naive", period=1, args=list(align.period=60))

# Traditional Estimate at 10 minute frequency 
rRealizedVariance(x=msft.real.cts[[1]], type="naive", period=10, args=list(align.period=60)) 

# Bartlett Kernel Estimate with minute aligned data at 20 lags 
rRealizedVariance(x=msft.real.cts[[1]], type="kernel", lags=20, args=list(align.period=60, type="Bartlett"))

# Cubic Kernel Estimate with second aligned data at 400 lags 
rRealizedVariance(x=msft.real.cts[[1]], type="kernel", lags=400, args=list(type="Cubic")) 

# Two-Timescale Estimate with minute aligned data at 10 subgrids 
rRealizedVariance(x=msft.real.cts[[1]], type="timescale", period=10, args=list(align.period=60))

# Subsample Average Estimate with second aligned data at 600 subgrids 
rRealizedVariance(x=msft.real.cts[[1]], type="avg", period=600) 


rc.zero(x=msft.real.cts[[1]], y=ge.real.cts[[1]], period=1)
rc.zero(x=msft.real.cts[[1]], y=ge.real.cts[[1]], period=60)

#java.graph("d:/dev/realizedDoc/Beta08/UsersManual/plots/Figure7A.jpg", format="JPEG", width=1400, height=1200)
#Figure 7A
test.zero <- rSignature(1:1200, x=msft.real.cts[[1]], y=ge.real.cts[[1]], type="zero", xscale=1/60)
plot(test.zero, ylab="% Zero", xlab="Sampling Frequency (Minutes)", main="MSFT | GE", sub=dates.example[[1]]) 
#dev.off()

#java.graph("d:/dev/realizedDoc/Beta08/UsersManual/plots/Figure7B.jpg", format="JPEG", width=1400, height=1200)
# Figure 7B
plot(rSignature((1:360)*5+1, x=msft.real.cts, y=ge.real.cts, xscale=1/60, iteration.funct="simpleIteration", iterations=1:5),
      ylab="Realized Covariance", xlab="Sampling Frequency (Minutes)", main="MSFT | GE", sub=paste(dates.example[[1]], dates.example[[5]], sep=" - "))
#dev.off()
      
#java.graph("d:/dev/realizedDoc/Beta08/UsersManual/plots/Figure8A.jpg", format="JPEG", width=1400, height=1200)
# Figure 8A
test.cov <- rSignature(1:1200,x=msft.real.cts[[1]], y=ge.real.cts[[1]], xscale=1/60) 
test.rect <- rSignature(1:600,msft.real.cts[[1]], ge.real.cts[[1]],type="kernel",args=list(type="rectangular"), xscale=1/30)
test.mth <- rSignature(1:600,msft.real.cts[[1]], ge.real.cts[[1]],type="kernel",args=list(type="mth"), xscale=1/30)
plot(test.cov, ylab="Realized Covariance", xlab="Minutes", main="GE | MSFT") 
lines(test.rect, col=3, lwd=1) 
lines(test.mth, col=4, lwd=2) 
axis(3, c(0,(1:5)*4), c("Lags:",as.character((1:5)*120)))
legend(13,.00015,c("Rectangular", "Mod TH"), lwd=c(1,2), col=c(3,4)) 

#dev.off()

#java.graph("d:/dev/realizedDoc/Beta08/UsersManual/plots/Figure8B.jpg", format="JPEG", width=1400, height=1200)
# Figure 8B
test.hy <- rSignature(1:600,msft.real.tts[[1]], ge.real.tts[[1]],type="hy",args=list(align.period=1), xscale=1/30) 
plot(test.cov, ylab="Realized Covariance", xlab="Minutes", main="GE | MSFT") 
lines(test.hy, col=2, lwd=2) 
axis(3, c(0,(1:5)*4), c("Tick Period:",as.character((1:5)*120))) 
legend(13,.00015,c("Hyashi-Yoshida"), lwd=c(2), col=c(2)) 

#dev.off()


# Traditional Estimate at highest frequency
rRealizedVariance(x=msft.real.cts[[1]], y=ge.real.cts[[1]], type="naive", period=1)

# Traditional Estimate at one minute frequency 
rRealizedVariance(x=msft.real.cts[[1]], y=ge.real.cts[[1]], type="naive", period=1, args=list(align.period=60)) 

# Traditional Estimate at 10 minute frequency 
rRealizedVariance(x=msft.real.cts[[1]], y=ge.real.cts[[1]], type="naive", period=10, args=list(align.period=60)) 

# Bartlett Kernel Estimate with minute aligned data at 20 lags 
rRealizedVariance(x=msft.real.cts[[1]], y=ge.real.cts[[1]], type="kernel", lags=20, args=list(align.period=60, type="Bartlett"))

# Cubic Kernel Estimate with second aligned data at 400 lags 
rRealizedVariance(x=msft.real.cts[[1]], y=ge.real.cts[[1]], type="kernel", lags=400, args=list(type="Cubic")) 

# Lead-Lag with one lag at one minute frequency 
rRealizedVariance(x=msft.real.cts[[1]], y=ge.real.cts[[1]], type="kernel", lags=1, args=list(align.period=60)) 
 
# Subsample Average Estimate with second aligned data at 600 subgrids 
rRealizedVariance(x=msft.real.cts[[1]], y=ge.real.cts[[1]], type="avg", period=600) 


# Traditional Estimate at highest frequency 
rRealizedVariance(x=merge(msft.real.cts[[1]], ge.real.cts[[1]]), type="naive", period=1)

# Traditional Estimate at 10 minute frequency 
rRealizedVariance(x=merge(msft.real.cts[[1]], ge.real.cts[[1]]), type="naive", period=10, args=list(align.period=60))

# Lead-Lag with one lag at one minute frequency> 
rRealizedVariance(x=merge(msft.real.cts[[1]], ge.real.cts[[1]]), type="kernel", lags=1, args=list(align.period=60))

# Subsample Average Estimate with second aligned data at 600 subgrids 
rRealizedVariance(x=merge(msft.real.cts[[1]], ge.real.cts[[1]]), type="avg", period=600)
 
 
 
 
# Traditional Estimate at highest frequency 
rRealizedVariance(x=msft.real.cts[[1]], y=ge.real.cts[[1]], type="naive", period=1, cor=T) 

# Traditional Estimate at one minute frequency 
rRealizedVariance(x=msft.real.cts[[1]], y=ge.real.cts[[1]], type="naive", period=1, args=list(align.period=60), cor=T)

# Traditional Estimate at 10 minute frequency 
rRealizedVariance(x=msft.real.cts[[1]], y=ge.real.cts[[1]], type="naive", period=10, args=list(align.period=60), cor=T) 

# Bartlett Kernel Estimate with minute aligned data at 20 lags 
rRealizedVariance(x=msft.real.cts[[1]], y=ge.real.cts[[1]], type="kernel", lags=20, args=list(align.period=60, type="Bartlett"), cor=T)

# Cubic Kernel Estimate with second aligned data at 400 lags 
rRealizedVariance(x=msft.real.cts[[1]], y=ge.real.cts[[1]], type="kernel", lags=400, args=list(type="Cubic"), cor=T) 

# Lead-Lag with one lag at one minute frequency 
rRealizedVariance(x=msft.real.cts[[1]], y=ge.real.cts[[1]], type="kernel", lags=1, args=list(align.period=60), cor=T) 

# Subsample Average Estimate with second aligned data at 600 subgrids 
rRealizedVariance(x=msft.real.cts[[1]], y=ge.real.cts[[1]], type="avg", period=600, cor=T) 


# Correlation Matrices
# Traditional Estimate at highest frequency 
rRealizedVariance(x=merge(msft.real.cts[[1]], ge.real.cts[[1]]), type="naive", period=1, cor=T)

# Traditional Estimate at 10 minute frequency
rRealizedVariance(x=merge(msft.real.cts[[1]], ge.real.cts[[1]]), type="naive", period=10, args=list(align.period=60), cor=T)

# Lead-Lag with one lag at one minute frequency 
rRealizedVariance(x=merge(msft.real.cts[[1]], ge.real.cts[[1]]), type="kernel", lags=1, args=list(align.period=60), cor=T) 

# Subsample Average Estimate with second aligned data at 600 subgrids >
rRealizedVariance(x=merge(msft.real.cts[[1]], ge.real.cts[[1]]), type="avg", period=600, cor=T) 




#java.graph("d:/dev/realizedDoc/Beta08/UsersManual/plots/Figure9A.jpg", format="JPEG", width=1400, height=1200)
# Figure 9A
cumm <- list() 
cumm[[1]] <- rCumSum(msft.real.cts[[1]], period=1, align.period=60) 
cumm[[2]] <- rCumSum(msft.real.cts[[1]], period=10, align.period=60) 
cumm[[3]] <- rCumSum(msft.real.cts[[1]], period=20, align.period=60) 
cumm[[4]] <- rCumSum(msft.real.cts[[1]], period=30, align.period=60) 
accum <- list() 
accum[[1]] <- rAccumulation(msft.real.cts[[1]], period=10, align.period=60) 
accum[[2]] <- rAccumulation(msft.real.cts[[1]], period=20, align.period=60) 
accum[[3]] <- rAccumulation(msft.real.cts[[1]], period=30, align.period=60)

par(mfrow=c(2,1)) 
plot(cumm[[1]], xlab="", ylab="Cumulative Ruturns", main="MSFT", sub=dates.example[[1]], type="p", col=16, lwd=2) 
lines(cumm[[2]], col=2, lwd=2) 
lines(cumm[[3]], col=3, lwd=2) 
lines(cumm[[4]], col=4, lwd=2) 
plot(accum[[1]], xlab="", ylab="Realized Accumulation", type="l",main="MSFT", sub=dates.example[[1]], col=2, lwd=2) 
lines(accum[[2]], col=3, lwd=2) 
lines(accum[[3]], col=4, lwd=2) 
#dev.off()


#java.graph("d:/dev/realizedDoc/Beta08/UsersManual/plots/Figure9B.jpg", format="JPEG", width=1400, height=1200)
# Figure 9B
par(mfrow=c(2,1))
plot(cumm[[2]], xlab="", ylab="Cumulative Ruturns", main="MSFT", sub=dates.example[[1]], type="p")
barplot(rMarginal(msft.real.cts[[1]], period=10, align.period=60)$y, main="Marginal Contribution Plot") 
#dev.off()



#java.graph("d:/dev/realizedDoc/Beta08/UsersManual/plots/Figure10A.jpg", format="JPEG", width=1400, height=1200)
# Figure 10A
cumm.ge <- list()
cumm.ge[[1]] <- rCumSum(ge.real.cts[[1]], period=1, align.period=60)
cumm.ge[[2]] <- rCumSum(ge.real.cts[[1]], period=10, align.period=60)

accum <- list() 
accum[[1]] <- rAccumulation(msft.real.cts[[1]], y=ge.real.cts[[1]], period=10, align.period=60)

par(mfrow=c(3,1))
plot(cumm[[1]], xlab="", ylab="Cumulative Ruturns", main="MSFT", sub=dates.example[[1]], type="p", col=16)
lines(cumm[[2]], col=2, lwd=3) 
plot(cumm.ge[[1]], xlab="", ylab="Cumulative Ruturns", main="GE", sub=dates.example[[1]], type="p", col=16) 
lines(cumm.ge[[2]], col=2, lwd=3)
plot(accum[[1]], xlab="", ylab="Realized Co-Accumulation", type="l",main="MSFT | GE", sub=dates.example[[1]], col=2) 
#dev.off()


#java.graph("d:/dev/realizedDoc/Beta08/UsersManual/plots/Figure10B.jpg", format="JPEG", width=1400, height=1200)
# Figure 10B
rScatterReturns(msft.real.cts[[1]],y=ge.real.cts[[1]], period=1, align.period=20,ylab="GE",xlab="MSFT",numbers=F) 
#dev.off()
par(mfrow=c(1,1))





# 
# Importing data without Finmetrics30 or Zivot, Yan
#
path <- "d:/dev/"
eur.usd.05.2007 <- read.table(paste(path, "eurusd.csv", sep=""), stringsAsFactors=F, sep=",")
usd.jpy.05.2007 <- read.table(paste(path, "usdjpy.csv", sep=""), stringsAsFactors=F, sep=",")
eur.jpy.05.2007 <- read.table(paste(path, "eurjpy.csv", sep=""), stringsAsFactors=F, sep=",")

eur.usd.05.2007[1:10,]  


getT <- function(x, dateStr,...) 
{ 
	y <- x[,3] 
	x[substring(y,1,10)==dateStr,] 
}

eur.usd.05.30.2007 <- getT(eur.usd.05.2007, "2007-05-30") 
usd.jpy.05.30.2007 <- getT(usd.jpy.05.2007, "2007-05-30") 
eur.jpy.05.30.2007 <- getT(eur.jpy.05.2007, "2007-05-30") 


midQuote <- function(x, bid.index = 4, ask.index = 5)
{
     (x[,bid.index] + x[,ask.index]	)/2
}

toMilliseconds <- function(x,...) { 
	ans <- 1000 * as.numeric(substring(x, 12,13)) * 60 * 60 + 
	       1000 * as.numeric(substring(x, 15,16)) * 60 + 
	       as.numeric(substring(x, 18,19)) * 1000 
	ans 
}



eur.usd.real <- realizedObject(list(data=midQuote(eur.usd.05.30.2007), 
                                    milliseconds=toMilliseconds(eur.usd.05.30.2007[,3])),
                                    makeReturns=T, cts=T, millisstart=0000, millisend=1000*24*60*60)
eur.jpy.real <- realizedObject(list(data=midQuote(eur.jpy.05.30.2007), 
                                    milliseconds=toMilliseconds(eur.jpy.05.30.2007[,3])),
                                    makeReturns=T, cts=T, millisstart=0000, millisend=1000*24*60*60)
usd.jpy.real <- realizedObject(list(data=midQuote(usd.jpy.05.30.2007), 
                                    milliseconds=toMilliseconds(usd.jpy.05.30.2007[,3])),
                                    makeReturns=T, cts=T, millisstart=0000, millisend=1000*24*60*60)


#java.graph("d:/dev/realizedDoc/Beta08/UsersManual/plots/Figure11A.jpg", format="JPEG", width=1400, height=1200)
# Figure 10A
test.sig <- rSignature(1:120, eur.usd.real, xscale=1/2, args=list(align.period=30)) 
test.bart <- rSignature(1:120, eur.usd.real, type="kernel",xscale=1/2, args=list(align.period=30, type="bartlett"))
plot(test.sig, ylab="Realized Variance", xlab="Minutes", main="Eur.Usd", sub="05/30/2007") 
lines(test.bart, col=2, lwd=2) 
axis(3, c(0,(1:10)*6), c("Lags:",as.character((1:10)*12))) 
#dev.off()


#java.graph("d:/dev/realizedDoc/Beta08/UsersManual/plots/Figure11B.jpg", format="JPEG", width=1400, height=1200)
# Figure 10B 
test.sig <- rSignature(1:60, x=eur.usd.real, y=eur.jpy.real, xscale=1/2, args=list(align.period=30), cor=T) 
test.bart <- rSignature(1:60, x=eur.usd.real, y=eur.jpy.real, type="kernel",xscale=1/2, args=list(align.period=30, type="bartlett"), cor=T) 
plot(test.sig, ylab="Realized Covariance", xlab="Minutes", main="Eur.Usd | Eur.Jpy", sub="05/30/2007") 
lines(test.bart, col=2, lwd=2) 
axis(3, c(0,(1:10)*3), c("Lags:",as.character((1:10)*6))) 
#dev.off()








