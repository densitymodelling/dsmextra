# Extrapolation tools for density surface models <img src="https://github.com/densitymodelling/dsmextra/blob/master/hex/dsmextra-hex.png?raw=true" height=200 align="right" class="logo"/>

<!-- badges: start -->
[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg?style=flat-square)](https://www.tidyverse.org/lifecycle/#maturing)
[![DOI](https://img.shields.io/badge/DOI-10.5281.zenodo.3529465-green.svg?style=flat-square)](https://doi.org/10.5281/zenodo.3529465)
<!-- badges: end -->

`dsmextra` provides a toolkit for quantifying and visualising extrapolation in ecological models (with a focus on density surface models, as implemented in package [dsm](https://cran.r-project.org/web/packages/dsm/index.html)) projected into novel environmental space. Currently, `dsmextra` defines extrapolation on the basis of two metrics: **(1) ExDet** (Mesgaran et al. 2014), and **(2) %N** (the percentage of data nearby, King & Zeng 2007). 

`dsmextra` offers a variety of numerical and graphical outputs, including summary plots and interactive maps created as [ggplot2](https://ggplot2.tidyverse.org/) and [html](https://rstudio.github.io/leaflet/) objects, respectively. Additional functionality (e.g. assessment methods for dynamic covariates) will be added in future releases.

The idea behind `dsmextra` is to aid ecologists, practitioners, and model end-users in identifying conditions (e.g. areas) under which predicted density surfaces may be prone to errors. In so doing, `dsmextra` may support:

+ Better-informed interpretations of (density surface) model outputs and their associated uncertainties.
+ Improvements to model development and covariate selection protocols.
+ Cost-effective allocation of future survey effort towards priority, data-poor areas.

### Getting started 

If you are just getting started with `dsmextra`, we recommend reading the tutorial [vignette](https://densitymodelling.github.io/dsmextra/articles/dsmextra-vignette.html), which provides a quick introduction to the package. You may also find the below technical report useful:

* Bouchet et al. (2019). [From here and now to there and then: Practical recommendations for extrapolating cetacean density surface models to novel conditions](https://research-repository.st-andrews.ac.uk/bitstream/handle/10023/18509/Denmod_ExtrapolationReport_final_Aug2019.pdf?sequence=1&isAllowed=y). CREEM technical report 2019-01, Centre for Research into Ecological & Environmental Modelling (CREEM), University of St Andrews, 59 p.

### Additional reading

* Mannocci et al. (2018). Assessing cetacean surveys throughout the mediterranean sea: A gap analysis in environmental space. *Scientific Reports* **8**, art3126. DOI: [10.1038/s41598-018-19842-9](https://www.nature.com/articles/s41598-018-19842-9).

* Mannocci et al. (2017). Extrapolating cetacean densities to quantitatively assess human impacts on populations in the high seas. *Conservation Biology* **31**, 601–614. DOI: [10.1111/cobi.12856](https://conbio.onlinelibrary.wiley.com/doi/full/10.1111/cobi.12856).

* Mesgaran et al. (2014). Here be dragons: A tool for quantifying novelty due to covariate range and correlation change when projecting species distribution models. *Diversity and Distributions* **20**, 1147–1159. DOI: [10.1111/ddi.12209](https://onlinelibrary.wiley.com/doi/full/10.1111/ddi.12209).

* Miller et al. (2013). Spatial models for distance sampling data: Recent developments and future directions. *Methods in Ecology and Evolution* **4**, 1001–1010. DOI: [10.1111/2041-210X.12105](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12105).

* King G & Zeng L (2007). When can history be our guide? The pitfalls of counterfactual inference. *International Studies Quarterly* **51**, 183–210. DOI: [10.1111/j.1468-2478.2007.00445.x](https://doi.org/10.1111/j.1468-2478.2007.00445.x).

### Acknowledgements

This R package was developed for the [DenMod project](https://synergy.st-andrews.ac.uk/denmod/) (Working group for the advancement of marine species density surface modelling), and was funded by OPNAV N45 and the SURTASS LFA Settlement Agreement, being managed by the U.S. Navy’s Living Marine Resources program under Contract No. N39430-17-C-1982.

### Installation

The latest development version can be installed from Github (requires the [remotes](https://github.com/r-lib/remotes) package):

```r
if (!require("remotes")) install.packages("remotes")
remotes::install_github("densitymodelling/dsmextra")
```

### Important note

**May 2020**

The [WhatIf](https://gking.harvard.edu/whatif) package, on which `dsmextra` depends, was removed from CRAN on 2020-03-07. Similarly, the [Zelig](https://zeligproject.org/) package (a dependency of WhatIf) was temporarily archived on 2020-03-24. As this may cause issues when first setting up `dsmextra`, a patch correction has been released under version 1.1.0 to ensure that `WhatIf` and `Zelig` are automatically pulled from Github during installation. Alternatively, both packages can be manually installed either from Github or from their respective archives (see [here](https://cran.r-project.org/web/packages/WhatIf/index.html) and [here](https://cran.r-project.org/web/packages/Zelig/index.html)). 

```r
if (!require("remotes")) install.packages("remotes")
remotes::install_github("IQSS/Zelig")
remotes::install_github("IQSS/WhatIf")
```

### Found a bug? Have a feature request?  

Please submit an issue or send a pull request to the [Github repository](https://github.com/densitymodelling/dsmextra/).
