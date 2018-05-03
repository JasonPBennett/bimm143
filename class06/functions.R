#' ---
#' title: "Class 6 Bioinformatics"
#' author: "Jason Patrick Bennett"
#' ---

# Functions

# Standard add function
add <- function(x, y=1) {
  # Sum x and y inputs
  x + y
}

# Rescaling function
rescale <- function(x) {
  rng <- range(x, na.rm = TRUE)
  rescaled <- (x - rng[1]) / (rng[2] - rng[1])
  return(rescaled)
}

# Second rescaling function
rescale2 <- function(x, na.rm = TRUE, plot = FALSE) {
  if(na.rm) {
    rng <- range(x, na.rm = na.rm)
  } else {
    rng <- range(x)
  }
  print("Hello")
  
  answer <- (x - rng[1]) / (rng[2] - rng[1])
  
  print("is it me you are looking for?")
  
  if(plot) {
    plot(answer, typ = "b", lwd = 4)
  }
  print("I can see it in...")
  
  return(answer)
}

# Test on a small example where you know the answer
rescale(1:10)

# Displays NA for all; Not optimal
rescale( c(1, 2, NA, 3, 10) )

# What should this do optimally?
rescale( c(1, 2, "string") )