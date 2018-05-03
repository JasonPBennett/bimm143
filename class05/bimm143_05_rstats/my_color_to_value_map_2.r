myColorMap2 <- function (x,
                         low.high = range(x),
                         palette = cm.colors(100)) {
  
  #Determine percent values of the high.low range
  percent <- ( (x - high.low[1])/(high.low[2] - high.low[1]) )
  
  #Find corresponding index position in the color palette
  #  Note: Catch for 0 but adding 1
  index <- round ( (length(palette) -1) * percent ) + 1
  
  return (palette[index])
}