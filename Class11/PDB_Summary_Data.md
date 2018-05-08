---
title: "Structurla Bioinformatics"
author: "Jason Patrick Bennett"
date: "May 8, 2018"
output:
  html_document:
    keep_md: true
---



## Structural Bioinformatics

Download the CSV file


```r
pdb.stats <- read.csv("Data Export Summary.csv")
```

Let's find the totals and their associated percentages


```r
percent <- (pdb.stats$Total / sum(pdb.stats$Total) ) * 100
names(percent) <- pdb.stats$Experimental.Method
percent
```

```
##               X-Ray                 NMR Electron Microscopy 
##         89.51253323          8.72181096          1.50770286 
##               Other        Multi Method 
##          0.17006317          0.08788979
```

