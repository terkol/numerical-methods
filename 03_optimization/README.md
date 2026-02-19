# Optimization

This folder contains small programs for unconstrained optimization and related eigenvalue-based “most unstable mode” searches. Topics include 1D bracketing minimization (golden-section search) and derivative-free multidimensional optimization (Nelder–Mead, Powell). One script (ballooning.py) uses eigenvalue computations to identify the mode with maximum growth rate, which is an optimization problem over discrete mode indices.

## Programs

### `golden_search.py`

Golden-section search for minimizing a unimodal function on an interval $[a,b]$.

The method maintains a bracket and evaluates $f$ at two interior points:

$r = (\sqrt{5}-1)/2$ (golden ratio complement)

$x_1 = a + (1-r)(b-a)$
$x_2 = a + r(b-a)$

and shrinks the interval depending on whether $f(x_1) < f(x_2)$, reusing one function evaluation each iteration. The interval length decreases by a constant factor per step.

Notes:

Correctness requires $f$ to be unimodal on $[a,b]$.

Currently the script only works for positive intervals. 

This script uses a polynomial–exponential objective involving an error-function approximation; in the current usage the search interval is $[0,2]$.

Dependencies: `numpy`, `matplotlib`.

Run: `python golden_search.py`

### `minimize.py`

Unconstrained 2D minimization using SciPy’s derivative-free optimizers.

Objective:
$f(x,y) = 10 e^{-(x-3)^2} + x^2 + xy + y^2 - 5 e^{-(y+2)^2}$.

The script runs two methods from different initial guesses:

- Nelder–Mead (simplex method)

- Powell (direction-set method)

These methods do not require gradients, but can converge to different local minima depending on initialization.

Dependencies: `numpy`, `scipy`.

Run: `python minimize.py`

### `ballooning.py`

Eigenvalue-based search for the most unstable mode in a discretized (ballooning-like) system.

Two comparisons are performed:

Uncoupled $2×2$ dispersion relation.
For each discrete $k_y$ mode, form a $2×2$ complex matrix and compute its eigenvalues $\omega$. The “growth rate” is typically taken as $\Im(\omega)$, and the script scans over modes to find the maximum growth rate:
$\max_{k_y} \Im(\omega(k_y))$.

Coupled $2M×2M$ sparse operator.
A block tridiagonal sparse matrix $L$ is assembled and the eigenmode with largest imaginary part is computed using ARPACK via SciPy:
$\omega, \psi : L\psi = \omega\psi$,
selecting the eigenvalue with maximal $\Im(\omega)$.

The script also visualizes:

- analytical vs numerical eigenvalues for the uncoupled case (real and imaginary parts),

- the spectral envelope of the selected mode,

- and its real-space representation via an inverse FFT constructed from a conjugate-symmetric spectrum.

Dependencies: `numpy`, `scipy`, `matplotlib`.

Run: `python ballooning.py`