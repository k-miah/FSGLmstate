# FSGLmstate <img src="https://raw.githubusercontent.com/k-miah/FSGLmstate/main/FSGLmstate.png" alt="Package Logo" align="right" height="250" />
[![](https://img.shields.io/badge/doi-10.48550/arXiv.2411.17394-yellow.svg)](https://doi.org/10.48550/arXiv.2411.17394)

`FSGLmstate` is an R package that performs variable selection via fused sparse-group lasso (FSGL) penalized multi-state models [(Miah et al., 2024)](https://doi.org/10.48550/arXiv.2411.17394).

## Abstract

In multi-state models based on high-dimensional data, effective modeling strategies are required to determine an optimal, ideally parsimonious model. 
In particular, linking covariate effects across transitions is needed to conduct joint variable selection. A useful technique to reduce model complexity is to address homogeneous covariate effects for distinct transitions. We integrate this approach to data-driven variable selection by extended regularization methods within multi-state model building. We propose the fused sparse-group lasso (FSGL) penalized Cox-type regression in the framework of multi-state models combining the penalization concepts of pairwise differences of covariate effects along with transition-wise grouping. For optimization, we adapt the alternating direction method of multipliers (ADMM) algorithm to Cox-type hazards regression in the multi-state setting. In a simulation study and application to acute myeloid leukemia (AML) data, we evaluate the algorithm's ability to select a sparse model incorporating relevant transition-specific effects and similar cross-transition effects. We investigate settings in which the combined penalty is beneficial compared to global lasso regularization.

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Features](#features)

## Installation

You can install the development package version `FSGLmstate` from GitHub with:

```R
# install.packages("devtools")
devtools::install_github("k-miah/FSGLmstate")
```

## Usage

Load the package in R with:

```R
library(FSGLmstate)
```

## Features

- `penalty_matrix_K()`: Generation of a penalty structure matrix for use in penalized regression incorporating lasso, fused and group-lasso penalties
- `fit.admm.fsgl.mstate()`: Alternating direction method of multipliers (ADMM) optimization for FSGL penalized multi-state models for fixed set of tuning parameters
- `gcv.fit.admm.fsgl.mstate()`: Alternating direction method of multipliers (ADMM) optimization for FSGL penalized multi-state models for optimal tuning parameters via generalized cross-validation (GCV) selection criterion

Bugs and issues can be reported at https://github.com/k-miah/FSGLmstate/issues.

## Contact

For any questions or feedback, please reach out to [k.miah@dkfz.de](mailto:k.miah@dkfz.de).


