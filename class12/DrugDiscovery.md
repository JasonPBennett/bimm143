---
title: "Drug Interactions"
author: "Jason Patrick Bennett"
date: "May 10, 2018"
output: 
  html_document:
    keep_md: true
---



# Drug Discovery

First load the Bio3D package and download our desired structure


```r
# Load the package
library(bio3d)

# Retrieve the HIV-1 Protease PDB file
file.name <- get.pdb("1hsg")
```

```
## Warning in get.pdb("1hsg"): ./1hsg.pdb exists. Skipping download
```

```r
# Save the file to the variable 'hiv'
hiv <- read.pdb(file.name)

# Quick look at the object
hiv
```

```
## 
##  Call:  read.pdb(file = file.name)
## 
##    Total Models#: 1
##      Total Atoms#: 1686,  XYZs#: 5058  Chains#: 2  (values: A B)
## 
##      Protein Atoms#: 1514  (residues/Calpha atoms#: 198)
##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
## 
##      Non-protein/nucleic Atoms#: 172  (residues: 128)
##      Non-protein/nucleic resid values: [ HOH (127), MK1 (1) ]
## 
##    Protein sequence:
##       PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYD
##       QILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNFPQITLWQRPLVTIKIGGQLKE
##       ALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTP
##       VNIIGRNLLTQIGCTLNF
## 
## + attr: atom, xyz, seqres, helix, sheet,
##         calpha, remark, call
```

**Question 1: What is the name of the two non-protein resid values?**
  Water (HOH) and MK1
  
**Question 2: What does resid correspond to and how would we get a listing of all residue values in this structure?**
  The position within the 1HSG where these things exist. You could get   a list of all the residue values by using: "hiv$atom$resid"
  
  
Now lets trim the file to make it easier to use for our data analysis


```r
# Select only the protein
prot <- trim.pdb(hiv, "protein")

# Select only the ligand (or the drug)
lig <- trim.pdb(hiv, "ligand")
```


And lets take these trimmed sections and create new pdb files from them


```r
# Create the protein file
write.pdb(prot, file = "1hsg_protein.pdb")

# Create the ligand file
write.pdb(lig, "1hsg_ligand.pdb")
```


Now use AutoDock to visualize the molecule

-- Attached Files --

Now lets inspect the docking results


```r
res <- read.pdb("all.pdbqt", multi = TRUE)

write.pdb(res, "results.pdb")
```


Pulling this up in VMD with the 1hsg.pdb file, we can see both the protein and our results superimposed on eachother. This allows us to *qualitatively* determine if our assumptions were correct.

---- Use VMD to compare the two ----

**Question 4: How good do the docks look qualitatively? Is the Crystal Binding Mode reproduced? Is it the best conformation according to AutoDock Vina?**
  They look good. ? ?

They look good, but lets determine *quantitatively* if they really are good.


```r
ori <- read.pdb("ligand.pdbqt")

rmsd(ori, res)
```

```
##  [1]  0.590 11.163 10.531  4.364 11.040  3.682  5.741  3.864  5.442 10.920
## [11]  4.318  6.249 11.084  8.929
```


**Question 5: How good are the docks quantitatively? Is the Crystal Binding Mode reproduced within 1 Angstrom RMSD for all atoms?**
  Our first dock is within 1 Angstrom, however the rest fall short.     Therefore, our first dock is the best.
  
**Question 6: How would you determine the RMSD for heavy atoms only (i.e. non Hydrogren atoms)?**
  See the code below.
  

```r
inds.res <- atom.select(res, string = "noh")

rmsd(lig, res$xyz[,inds.res$xyz])
```

```
##  [1]  7.073  8.765  9.086  7.101  9.798  6.967  7.556  7.219  7.362  9.201
## [11]  7.295  7.213 10.257 10.840
```

