Homework
================
Xuqian Tan
1/30/2019

Original Code
-------------

``` r
library(bio3d)
s1 <- read.pdb("4AKE") # kinase with drug
```

    ##   Note: Accessing on-line PDB file

``` r
s2 <- read.pdb("1AKE") # kinase no drug
```

    ##   Note: Accessing on-line PDB file
    ##    PDB has ALT records, taking A only, rm.alt=TRUE

``` r
s3 <- read.pdb("1E4Y") # kinase with drug
```

    ##   Note: Accessing on-line PDB file

``` r
s1.chainA <- trim.pdb(s1, chain="A", elety="CA")
s2.chainA <- trim.pdb(s2, chain="A", elety="CA")
s3.chainA <- trim.pdb(s3, chain="A", elety="CA")

s1.b <- s1.chainA$atom$b
s2.b <- s2.chainA$atom$b
s3.b <- s3.chainA$atom$b

plotb3(s1.b, sse=s1.chainA, typ="l", ylab="Bfactor")
```

![](BIMM143_Homework_Xuqian_Tan_files/figure-markdown_github/unnamed-chunk-1-1.png)

``` r
plotb3(s2.b, sse=s2.chainA, typ="l", ylab="Bfactor")
```

![](BIMM143_Homework_Xuqian_Tan_files/figure-markdown_github/unnamed-chunk-1-2.png)

``` r
plotb3(s3.b, sse=s3.chainA, typ="l", ylab="Bfactor")
```

![](BIMM143_Homework_Xuqian_Tan_files/figure-markdown_github/unnamed-chunk-1-3.png)

My function & Documentation:
----------------------------

#### Function name: Analysis

#### Description: This function takes in a stirng of protein name, analyzes protein drug interactions by reading in any protein PDB data and outputs a plot for the specified protein

#### Input: name -&gt; a string of the name of the interested protein

#### Output: a plot of protein drug interaction of the interested protein

``` r
Analysis <- function(name) {
  s <- read.pdb(name)
  s.chainA <- trim.pdb(s, chain="A", elety="CA")
  s.b <- s.chainA$atom$b
  plotb3(s.b, sse=s.chainA, typ="l", ylab="Bfactor")
}
```

### test

``` r
names <- c("4AKE", "1AKE", "1E4Y")
for (i in names){
  Analysis(i)
}
```

    ##   Note: Accessing on-line PDB file

    ## Warning in get.pdb(file, path = tempdir(), verbose = FALSE): /var/folders/
    ## m7/klc0vqy504s_v0498wvlm09r0000gn/T//RtmpVsrF75/4AKE.pdb exists. Skipping
    ## download

    ##   Note: Accessing on-line PDB file

    ## Warning in get.pdb(file, path = tempdir(), verbose = FALSE): /var/folders/
    ## m7/klc0vqy504s_v0498wvlm09r0000gn/T//RtmpVsrF75/1AKE.pdb exists. Skipping
    ## download

![](BIMM143_Homework_Xuqian_Tan_files/figure-markdown_github/unnamed-chunk-3-1.png)

    ##    PDB has ALT records, taking A only, rm.alt=TRUE

    ##   Note: Accessing on-line PDB file

    ## Warning in get.pdb(file, path = tempdir(), verbose = FALSE): /var/folders/
    ## m7/klc0vqy504s_v0498wvlm09r0000gn/T//RtmpVsrF75/1E4Y.pdb exists. Skipping
    ## download

![](BIMM143_Homework_Xuqian_Tan_files/figure-markdown_github/unnamed-chunk-3-2.png)![](BIMM143_Homework_Xuqian_Tan_files/figure-markdown_github/unnamed-chunk-3-3.png)
