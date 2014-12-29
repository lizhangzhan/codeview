library(latticeExtra)
library(KFAS)



xyplot(EuStockMarkets[,1]) +
  layer(panel.tskernel(x, y, width = 100, col = 1, sides = 1)) +
  layer(panel.tskernel(x, y, width = 20, col = 2, sides = 1))

x = EuStockMarkets[,1]

plot(x)






library(KFAS)
library(tseries)
library(timeSeries)
library(zoo)
library(quantmod)

getDailyPrices = function( tickerSym, startDate, endDate )
{
  prices = get.hist.quote( instrument = tickerSym, start = startDate, end = endDate,
                           quote="AdjClose", provider="yahoo", 
                           compression="d",  quiet=T)
  
  prices.ts = ts(prices)
  return( prices.ts )
}

kalmanFilter = function( x )
{
  t = x
  if (class(t) != "ts") {
    t = ts(t)
  }
  ssModel = structSSM( y = t, distribution="Gaussian")
  ssFit = fitSSM(inits=c(0.5*log(var(t)), 0.5*log(var(t))), model = ssModel )
  kfs = KFS( ssFit$model, smoothing="state", nsim=length(t))
  vals = kfs$a
  lastVal = vals[ length(vals)]
  return(lastVal)
}

Start = "2011-01-01"
End   = "2012-12-31"
SandP = "^GSPC"

windowWidth = 20
tsLength = 100

SAndP.ts = getDailyPrices( SandP, Start, End )
SAndP.ts = SAndP.ts[1:tsLength]
SAndP.smoothed = rollapply( data=SAndP.ts, width=windowWidth, FUN=kalmanFilter)

par(mfrow=c(1,1))
prices = coredata( SAndP.ts[windowWidth:length(SAndP.ts)])
plot(prices, col="blue", type="l")
lines(coredata(SAndP.smoothed), col="magenta")
par(mfrow=c(1,1))





library(mFilter)
data(unemp)
opar <- par(no.readonly=TRUE)
unemp.bk <- bkfilter(unemp)
plot(unemp.bk)
unemp.bk1 <- bkfilter(unemp, drift=TRUE)
unemp.bk2 <- bkfilter(unemp, pl=8,pu=40,drift=TRUE)
unemp.bk3 <- bkfilter(unemp, pl=2,pu=60,drift=TRUE)
unemp.bk4 <- bkfilter(unemp, pl=2,pu=40,drift=TRUE)
par(mfrow=c(2,1),mar=c(3,3,2,1),cex=.8)
plot(unemp.bk1$x,
     main="Baxter-King filter of unemployment: Trend, drift=TRUE",
     col=1, ylab="")
lines(unemp.bk1$trend,col=2)
lines(unemp.bk2$trend,col=3)
lines(unemp.bk3$trend,col=4)
lines(unemp.bk4$trend,col=5)
legend("topleft",legend=c("series", "pl=2, pu=32", "pl=8, pu=40",
                          "pl=2, pu=60", "pl=2, pu=40"), col=1:5, lty=rep(1,5), ncol=1)


## library(mFilter)
data(unemp)
opar <- par(no.readonly=TRUE)
unemp.bw <- bwfilter(unemp)
plot(unemp.bw)
unemp.bw1 <- bwfilter(unemp, drift=TRUE)
unemp.bw2 <- bwfilter(unemp, freq=8,drift=TRUE)
unemp.bw3 <- bwfilter(unemp, freq=10, nfix=3, drift=TRUE)
unemp.bw4 <- bwfilter(unemp, freq=10, nfix=4, drift=TRUE)
par(mfrow=c(2,1),mar=c(3,3,2,1),cex=.8)
plot(unemp.bw1$x,
     main="Butterworth filter of unemployment: Trend,
drift=TRUE",col=1, ylab="")
lines(unemp.bw1$trend,col=2)
lines(unemp.bw2$trend,col=3)
lines(unemp.bw3$trend,col=4)
lines(unemp.bw4$trend,col=5)
legend("topleft",legend=c("series", "freq=10, nfix=2",
                          "freq=8, nfix=2", "freq=10, nfix=3", "freq=10, nfix=4"),
       col=1:5, lty=rep(1,5), ncol=1)





data(unemp)
opar <- par(no.readonly=TRUE)
unemp.cf <- cffilter(unemp)
plot(unemp.cf)
unemp.cf1 <- cffilter(unemp, drift=TRUE, root=TRUE)
unemp.cf2 <- cffilter(unemp, pl=8,pu=40,drift=TRUE, root=TRUE)
unemp.cf3 <- cffilter(unemp, pl=2,pu=60,drift=TRUE, root=TRUE)
unemp.cf4 <- cffilter(unemp, pl=2,pu=40,drift=TRUE, root=TRUE,theta=c(.1,.4))
par(mfrow=c(2,1),mar=c(3,3,2,1),cex=.8)
plot(unemp.cf1$x,
     main="Christiano-Fitzgerald filter of unemployment: Trend \n root=TRUE,drift=TRUE",
     col=1, ylab="")
lines(unemp.cf1$trend,col=2)
lines(unemp.cf2$trend,col=3)
lines(unemp.cf3$trend,col=4)
lines(unemp.cf4$trend,col=5)
legend("topleft",legend=c("series", "pl=2, pu=32", "pl=8, pu=40", "pl=2, pu=60",
                          "pl=2, pu=40, theta=.1,.4"), col=1:5, lty=rep(1,5), ncol=1)



