# pasta: Pattern Analysis for Spatial Omics daTA

This repository contains code for the `pasta` vignette accompanying the [pasta preprint](
https://doi.org/10.48550/arXiv.2412.01561).

## Disclaimer

This vignette should be a dynamic resource for people interested in spatial statistics for spatial omics data. The content is therefore subject to change. If you wish to contribute to this resource, feel free to open an issue or contact us via email.

## Installing all necessary R packages

`spatialFDA` has to be installed via `devtools::install_github("mjemons/spatialFDA")` all other R packages can be installed using `renv::restore()`.

## Using `renv` and Python

It can happen that `renv` creates a virtual environment that does not locate the correct Python version and it will not be able to find the shared libraries.

To get a (hopefully) working environment: 

- Select your standart Python installation with `renv::use_python()`. It has to correspond to the Python version in the lock file (`3.10.12`).
- This will create a virtual environment with the Python binary pointing to the selected version.
- You can then use `renv::restore` to install both R and Python packages in the `renv` (virtual) environment.

Disclaimer: This was so far only tested on `3.10.12 | packaged by conda-forge` and `macOS Sonoma 14.7.1 - aarch64-apple-darwin20 (64-bit)`.

