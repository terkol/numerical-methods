# Numerical methods and scientific computing mini-programs

This repository contains a collection of small programs from a course in numerical methods / scientific computing. The code is organized by topic and includes short demonstrations of numerical accuracy, stability, and standard algorithms used in scientific computing.

Languages: Python and Fortran.

## Folder overview

### `01_floating_point_and_errors`

Floating-point effects and numerical error: catastrophic cancellation, rounding-error accumulation in summation, and truncation vs rounding tradeoffs in numerical differentiation.

### `02_root_finding_and_nonlinear_systems`

Scalar root finding (bisection, Newton), nonlinear systems (Broyden), and related examples.

### `03_optimization`

1D bracketing minimization (golden-section search), derivative-free multivariate minimization (SciPy), and an eigenvalue-based “most unstable mode” search.

### `04_linear_algebra_and_eigenproblems`

Dense linear solves (LU via LAPACK), least-squares fitting, eigenvalue algorithms, and polynomial roots via a companion-matrix eigenvalue formulation.

### `05_interpolation_and_splines`

Polynomial interpolation (Lagrange/Newton forms) and spline methods (B-splines, cubic splines).

### `06_numerical_integration`

Quadrature methods including Romberg integration (Richardson extrapolation of trapezoids), Simpson refinement, and Gauss–Laguerre quadrature.

### `07_odes_and_boundary_value_problems`

Ordinary differential equations (ODEs): RK4 for initial value problems and shooting for boundary value problems.

### `08_fourier_and_spectral_methods`

FFT-based signal analysis (frequency estimation, spectrograms) and spectral time stepping via Fourier transforms.

### `09_probability_and_statistics`

Monte Carlo and statistical demonstrations: bootstrap confidence intervals, central limit theorem behavior (including heavy-tailed counterexamples), and chi-squared tests for randomness.

### `10_parameter_estimation_and_fitting`

Parameter estimation from data: MCMC likelihood fitting and nonlinear least-squares fitting of a Morse potential.

### `11_quantum_applications`

A finite-difference Schrödinger eigenvalue problem on a 1D interval (particle in a box and a simple perturbation).

Each topic folder contains a `README.md` describing the programs in that folder and how to run them.

## Running the code

### Python

Most Python scripts require `numpy`. Some additionally require `matplotlib`, `scipy`, or `sympy`. A minimal installation that covers most scripts is:

`pip install numpy matplotlib scipy sympy`

One script uses optional packages (e.g. `pyfftw` for FFTW-backed FFTs). Topic folder READMEs list any extra dependencies.

Run a script directly, for example:

`python src/02_root_finding_and_nonlinear_systems/broydens_method.py`

### Fortran

Standalone Fortran programs can be compiled with:

`gfortran file.f90 -O2 -o prog && ./prog`

Some Fortran programs call LAPACK/BLAS routines (e.g. `DGESV`, `DGEEV`) and must be linked accordingly:

`gfortran file.f90 -O2 -llapack -lblas -o prog && ./prog`

(Exact link flags can vary across systems depending on how BLAS/LAPACK are installed.)

## Numerical stability notes

Several programs are designed to illustrate numerical failure modes (for example underflow in product likelihoods or catastrophic cancellation). In those cases the relevant folder README explains what is being demonstrated and what a numerically stable alternative would be.

## Data files

Some scripts load data files such as `sample*.dat`, `linear_data*.txt`, or `morse_data.txt`. These files are stored in the data folder. The corresponding folder README indicates the expected input format.
