# FSGLmstate

`FSGLmstate` is an R package that performs variable selection for fused sparse-group lasso (FSGL) penalized multi-state models (Miah et al., 2024).

## Abstract

In multi-state models based on high-dimensional data, effective modeling strategies are required to determine an optimal, ideally parsimonious model. 
In particular, linking covariate effects across transitions is needed to conduct joint variable selection. A useful technique to reduce model complexity is to address homogeneous covariate effects for distinct transitions. We integrate this approach to data-driven variable selection by extended regularization methods within multi-state model building. We propose the fused sparse-group lasso (FSGL) penalized Cox-type regression in the framework of multi-state models combining the penalization concepts of pairwise differences of covariate effects along with transition grouping. For optimization, we adapt the alternating direction method of multipliers (ADMM) algorithm to transition-specific hazards regression in the multi-state setting. In a simulation study and application to acute myeloid leukemia (AML) data, we evaluate the algorithm's ability to select a sparse model incorporating relevant transition-specific effects and similar cross-transition effects. We investigate settings in which the combined penalty is beneficial compared to only fused or grouped regularization.

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Features](#features)

## Installation

You can install the package FSGLmstate from GitHub with:

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

- Multi-state partial log-likelihood function with first and second derivatives
- Multi-state Cox estimation algorithms gradient ascent & Newton-Raphson based on long format data
- Alternating direction method of multipliers (ADMM) for FSGL penalized multi-state models
- Choice of optimal tuning parameter by generalized cross-validation (GCV)

## Contact

For any questions or feedback, please reach out to [k.miah@dkfz.de](mailto:k.miah@dkfz.de).


