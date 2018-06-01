---
title: "Bioinformatics Class 7"
author: "Jason Patrick Bennett"
date: "April 24, 2018"
output:
  html_document:
    keep_md: yes
---



## Functions: *Part Two*

Here I am going to revsit our function from class 6: **rescale**.


```r
source("http://tinyurl.com/rescale-R")
```

Lets see if we can use this function!


```r
rescale(1:10)
```

```
##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
##  [8] 0.7777778 0.8888889 1.0000000
```

Looks good! Lets break it!


```r
#rescale( c(1:10, "break") )
```

Lets try the new **rescale2()** function.


```r
#rescale2( c(1:10, "brokenButBetter") )
```

## Write a NA checking function

Here we write a new function to check for NAs in two inputs


```r
x <- c( 1, 2, NA, 3, NA)
y <- c(NA, 3, NA, 3, 4)
```



```r
is.na(x)
```

```
## [1] FALSE FALSE  TRUE FALSE  TRUE
```



```r
which (is.na(x) )
```

```
## [1] 3 5
```



```r
sum ( is.na(x) )
```

```
## [1] 2
```



```r
sum ( is.na(x) & is.na(y) )
```

```
## [1] 1
```



```r
which ( is.na(x) & is.na(y) )
```

```
## [1] 3
```



```r
both_na(x, y)
```

```
## [1] 1
```



```r
x <- c(NA, NA, NA)
y1 <- c(1, NA, NA)
y2 <- c(1, NA, NA, NA)

#both_na(x, y2)
```



```r
#both_na2(x, y2)
```



```r
both_na3(x, y1)
```

```
## Found 2 NA's at position(s):2, 3
```

```
## $number
## [1] 2
## 
## $which
## [1] 2 3
```


## Another function example: gene intersection


```r
df1
```

```
##     IDs exp
## 1 gene1   2
## 2 gene2   1
## 3 gene3   1
```

```r
df2
```

```
##     IDs exp
## 1 gene2  -2
## 2 gene4  NA
## 3 gene3   1
## 4 gene5   2
```

```r
x <- df1$IDs
y <- df2$IDs
```



```r
x
```

```
## [1] "gene1" "gene2" "gene3"
```

```r
y
```

```
## [1] "gene2" "gene4" "gene3" "gene5"
```



```r
#intersect(x, y)

x %in% y
```

```
## [1] FALSE  TRUE  TRUE
```



```r
y %in% x
```

```
## [1]  TRUE FALSE  TRUE FALSE
```


Now we can access the genes we want with these indices!



```r
x[x %in% y]
```

```
## [1] "gene2" "gene3"
```

```r
y[y %in% x]
```

```
## [1] "gene2" "gene3"
```


Can make these columns of the same object with **cbind()**



```r
cbind( x[ x %in% y], y[ y %in% x] )
```

```
##      [,1]    [,2]   
## [1,] "gene2" "gene2"
## [2,] "gene3" "gene3"
```



```r
gene_intersect(x, y)
```

```
##      [,1]    [,2]   
## [1,] "gene2" "gene2"
## [2,] "gene3" "gene3"
```



```r
gene_intersect2(df1, df2)
```

```
##     IDs exp df2[df2$IDs %in% df1$IDs, "exp"]
## 2 gene2   1                               -2
## 3 gene3   1                                1
```


Lets try the **merge()** function for this job



```r
merge(df1, df2, by="IDs")
```

```
##     IDs exp.x exp.y
## 1 gene2     1    -2
## 2 gene3     1     1
```

