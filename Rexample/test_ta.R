library(quantmod)
library(xts)
start = '2010-01-01'
end = '2014-03-24'
ticker = 'AAPL'
f = getSymbols(ticker, src = 'yahoo', from = start, to = end, auto.assign=F)

require(TTR)

chartSeries(
  as.vector(f[,'AAPL.Close']),
  theme = chartTheme("white"),
  TA = c(addSMA(20),addSMA(50))
)

a <- f[,'AAPL.Close']

SMA(as.vector(a))


chartSeries(
  f[,'AAPL.Close'],
  theme = chartTheme("white"),
  TA = c(addSMA(20),addSMA(50),addBBands(),addRSI())
)

'20101020 00:00:00'



require(xts)
require(magrittr)

data(sample_matrix)
class(sample_matrix)
str(sample_matrix)
as.data.frame(sample_matrix) %>% rownames %>% head
