map.colors <- function (value,high.low,palette) {
  proportion <- ((value-high.low[1])/(high.low[2]-high.low[1]))
  index <- round ((length(palette)-1)*proportion)+1
  return (palette[index])
}

#High.low is a vector of two values that the user needs to create: the highest value and the lowest value
# in the expression column