
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

[![minimal R
version](https://img.shields.io/badge/R%3E%3D-4.3.0-brightgreen.svg)](https://cran.r-project.org/)
[![Licence](https://img.shields.io/github/license/mashape/apistatus.svg)](http://choosealicense.com/licenses/mit/)
[![R build
status](https://github.com/OJWatson/hrpup/workflows/R-CMD-check/badge.svg)](https://github.com/OJWatson/hrpup/actions)
[![CodeFactor](https://www.codefactor.io/repository/github/OJWatson/hrpup/badge)](https://www.codefactor.io/repository/github/OJWatson/hrpup)
<!-- badges: end -->

## Research compendium for hrp2 risk mapping

This is a working R compendium (think R package but for reproducible
analysis). A good overview on research compendiums, see the [R for
Reproducible Research](https://annakrystalli.me/rrresearch/index.html)
course.

### Installation

    git clone https://github.com/OJWatson/hrpup.git
    cd hrpup
    open hrpup.Rproj

Next, if `renv` has been used in this repository (look out for
`renv.lock`) then use `renv::restore` to set up package dependencies.
Otherwise `devtools::install_dev_deps()` will install all required
packages, as specified in the Imports in DESCRIPTION.

### Overview

The structure within analysis is as follows:

    R/                            # Packaged R functions 

    analysis/
        |
        ├── 01_xxxxx /            # analysis scripts used for generating figures
        |
        ├── plots/                # location of figures produced by the analysis scripts
        |
        ├── tables/               # location of tables produced by the analysis scripts
        |
        ├── data_raw/             # data obtained from elsewhere and treated read-only    
        |
        ├── data_derived/         # intermediate data generated during the analysis
        |
        ├── data_out/             # data outputs produced for external partners

Any analysis scripts with “X\_” in the name are used to format raw data
shared with us for this project that could not be included in the
repository and only summaries of the data as needed for the modelling
could be included. These scripts are still included to show how the raw
data was processed and the summaries saved in `data_raw` for use later
on. In this way, we ensure transparency and reproducibility.

Analysis scripts are to be run in the numbered order they are included.
If there are shared numbers, then any order of those scripts works.

### Compendium DOI:

<https://doi.org/X/X>

The files at the URL above will generate the results as found in the
publication.

### The R package

This repository is organized as an R package. There are only a few R
functions exported in this package - the majority of the R code is in
the analysis directory. The R package structure is here to help manage
dependencies, to take advantage of continuous integration, and so we can
keep file and data management simple. For any R packages that are used
frequently in this repository, they are documented in `R/` and are used
in the analysis folder using `devtools::load_all()`.

To download the package source as you see it on GitHub, for offline
browsing, use this line at the shell prompt (assuming you have Git
installed on your computer):

``` r
git clone https://github.com/OJWatson/hrpup.git
```

Once the download is complete, open the `hrpup.Rproj` in RStudio to
begin working with the package and compendium files. We will endeavour
to keep all package dependencies required listed in the DESCRIPTION.

<!-- To add this once all the analysis is done -->

In addition, once analysis is completed, we will use `renv` to track
package dependencies for reproducibility. Please use `renv::restore` to
restore the state of the project and see
<https://rstudio.github.io/renv/articles/renv.html> for more
information.

### Licenses

Code: [MIT](http://opensource.org/licenses/MIT) year: 2023, copyright
holder: OJ Watson

Data: [CC-0](http://creativecommons.org/publicdomain/zero/1.0/)
attribution requested in reuse
