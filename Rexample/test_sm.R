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






