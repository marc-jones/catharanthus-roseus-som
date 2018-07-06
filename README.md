# **Catharanthus roseus** self-organizing map clustering

The files in this repository recreate figure 2 from the paper "An NPF transporter exports a central monoterpene indole alkaloid intermediate from the vacuole"

The transcriptome data is the same as that available from the [Medicinal Plant Genomics](http://medicinalplantgenomics.msu.edu/pub/data/MPGR/Catharanthus_roseus) website, formatted such that R can more easily read it in.

To run the script, set the directory to the working directory within R and simply run:

```R
source(som.R)
```

This is the recommended way of running the script as it will prompt the user to input the number of clusters to select for further analysis.

The script requires the `kohonen` and `fBasics` packages. Specifically, it requires version 2.0.19 of the `kohonen` package, which can be installed as follows:

```R
install.packages('devtools') # If devtools is not already installed
library(devtools)
install_version('kohonen', version='2.0.19')
```

It is reasonably trivial to update this script to use the newer versions of the `kohonen` package, but for the perposes of reproducability, the older version is used here.

# References

[Payne, R. M. E. et al. An NPF transporter exports a central monoterpene indole alkaloid intermediate from the vacuole. Nat. Plants 3, 16208 (2017).](https://doi.org/10.1038/nplants.2016.208)
