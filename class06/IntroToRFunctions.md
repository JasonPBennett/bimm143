---
title: "Introduction to R Functions"
author: "Jason Patrick Bennett"
date: "June 1, 2018"
output:
  html_document:
    keep_md: true
---



## R Functions

Lets first try a standard add function:


```r
# Standard add function
add <- function(x, y=1) {
  # Sum x and y inputs
  x + y
}
```

And now, a function to rescale data:


```r
# Rescaling function
rescale <- function(x) {
  rng <- range(x)#, na.rm = TRUE)
  rescaled <- (x - rng[1]) / (rng[2] - rng[1])
  return(rescaled)
}
```

Now we'll try a more advanced rescale function that also includes a plot depending on a logical parameter passed to the function:


```r
# Second rescaling function
rescale2 <- function(x, na.rm = TRUE, plot = FALSE) {
  if(na.rm) {
    rng <- range(x, na.rm = TRUE)
  } else {
    rng <- range(x)
  }
  
  answer <- (x - rng[1]) / (rng[2] - rng[1])
  
  if(plot) {
    plot(answer, typ = "b", lwd = 4)
  }

  return(answer)
}
```

Lets test the first function now:


```r
# Test on a small example where you know the answer
rescale(1:10)
```

```
##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
##  [8] 0.7777778 0.8888889 1.0000000
```


```r
# Displays NA for all; Not optimal
rescale( c(1, 2, NA, 3, 10) )
```

```
## [1] NA NA NA NA NA
```

This second rescale result was not exactly what we were looking for. What can we do to improve this?

  - Remove NA from the dataset before trying to rescale! (Uncomment the na.rm = TRUE      segment of the rescale function).


```r
rescale2( c(1, 2, NA, 3, 10))
```

```
## [1] 0.0000000 0.1111111        NA 0.2222222 1.0000000
```
