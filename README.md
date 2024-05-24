# FSGLmstate

`FSGLmstate` is an R package that performs variable selection for fused sparse-group lasso (FSGL) penalized multi-state models (Miah et al., 2024).

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

- Multi-state partial log-likelihood function with derivatives
- Multi-state Cox estimation algorithms: Gradient ascent & Newton-Raphson based on long format data
- Alternatin direction method of multipliers (ADMM) for FSGL penalized multi-state models
- Choice of optimal tuning parameter by generalized cross-validation (GCV)

## Contact

For any questions or feedback, please reach out to [k.miah@dkfz.de](mailto:k.miah@dkfz.de).


