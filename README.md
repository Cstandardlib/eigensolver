# eigensolver: A C++/Eigen Implementation of the LOBPCG Algorithm

Eigensolver is a C++ library that provides an implementation of the Locally Optimal Block Preconditioned Conjugate Gradient (LOBPCG) algorithm for solving large-scale sparse eigenvalue problems. The library is built on top of the Eigen framework, which is a high-level C++ library for linear algebra.

## Table of Contents

- [eigensolver: A C++/Eigen Implementation of the LOBPCG Algorithm](#eigensolver-a-ceigen-implementation-of-the-lobpcg-algorithm)
  - [Table of Contents](#table-of-contents)
  - [Introduction](#introduction)
  - [Features](#features)
  - [Dependencies](#dependencies)
  - [Installation](#installation)
  - [Algorithm](#algorithm)
    - [Examples](#examples)
  - [Contributing](#contributing)

## Introduction

Eigensolver is designed to find a few eigenvalues and eigenvectors of a large sparse matrix, which is a common task in scientific computing and engineering. The LOBPCG method is a powerful numerical algorithm for this purpose, and we implement a version to seek smallest eigenpairs.

## Features

- Solves eigenvalue problems for large sparse matrices. Supports both standard and general symmetric eigenvalue problems. 
- Utilizes the power of the Eigen library for linear algebra operations.
- Provides a standalone implementation and an easy-to-use API for integrating into existing C++ projects, based on Eigen data structures and input matrices organized as [mtx](https://math.nist.gov/MatrixMarket/) files.
- Use [xmake](https://xmake.io) as build tools.

## Dependencies

- **Eigen**: A high-level C++ library for linear algebra.
- **fast_matrix_market**: A library for reading and writing Matrix Market files. Provides interface for eigen support.

## Installation

Eigensolver requires a C++ compiler with support for C++11 or higher and the Eigen library. 

The Eigen library is included in the source code, and fast_matrix_market is organized as a submodule.


To install and set up Eigensolver, follow these steps:

1.  Clone the Eigensolver repository:
```bash
git clone --recurse-submodules https://github.com/Cstandardlib/eigensolver.git
```

2. Build your projects with [xmake](https://xmake.io).
```bash
$ xmake build main
$ xmake run main
```

<!-- ## Usage

To use Eigensolver in your project, simply include the relevant headers and follow the provided API:

```cpp
#include "eigensolver/lobpcg.h"

int main() {
 // Define your sparse matrix A and initial guess X
 // ...

 // Set up the LOBPCG solver
 Eigensolver::LOBPCG solver;
 solver.setInputMatrix(A);
 solver.setInitialGuess(X);

 // Run the solver
 solver.solve();

 // Retrieve the eigenvalues and eigenvectors
 auto eigenvalues = solver.getEigenvalues();
 auto eigenvectors = solver.getEigenvectors();

 return 0;
}
``` -->

## Algorithm
The LOBPCG algorithm is an iterative method used to find selected eigenvalues and eigenvectors of a large sparse matrix. It is particularly effective when only a small number of eigenvalues and eigenvectors are required. The algorithm can be summarized in the following steps:

1. Start with an initial guess for the eigenvectors.
2. Perform a block conjugate gradient iteration to find a more accurate approximation.
3. Apply a preconditioner if provided.
4. Compute the Rayleigh quotient to estimate the eigenvalues.
5. Repeat the process until convergence.
### Examples
Eigensolver comes with several examples that demonstrate how to use the library to solve different types of eigenvalue problems. You can find these examples in the examples/ directory of the repository.

## Contributing
Contributions to Eigensolver are welcome! If you find a bug, have an idea for a new feature, or want to improve the documentation, please open an issue or submit a pull request.

<!-- ## License
Eigensolver is released under the MIT License. Feel free to use it in your commercial or personal projects. -->

