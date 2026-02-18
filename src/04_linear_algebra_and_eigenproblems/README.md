# Linear algebra and eigenproblems

This folder contains small programs illustrating core numerical linear algebra tasks: solving linear systems, least-squares fitting, eigenvalue computation, and the sensitivity of eigenvalues to perturbations.

## Programs

### `lu_solve_residual.f90`

Solves a dense linear system

$A x = b$

using LAPACK DGESV (LU factorization with partial pivoting), then computes the residual

$r = A x - b$

and reports a residual norm $|r|$.

Notes:

DGESV overwrites its inputs, so the code makes copies of $A$ and $b$ before calling LAPACK.

The residual is computed using the original (unmodified) $A$ and $b$.

Build/run (requires LAPACK/BLAS):
`gfortran lu_solve_residual.f90 -O2 -llapack -lblas -o lu_solve_residual && ./lu_solve_residual < input.txt`

Input format (via stdin):

first line: $N$

next $N$ lines: rows of $A$

final line: entries of $b$

### `noisy_matrix.py`

Linear least-squares parameter estimation for a model of the form

$y(t) \approx \alpha + \beta t + \gamma^2 \cos(t)$.

The script constructs the design matrix $A$ with rows

$[1,; t,; \cos(t)]$

and estimates parameters by solving the normal equations

$(A^T A),x = A^T y$,

and also via the Moore–Penrose pseudoinverse

$x = A^+ y$.

It then estimates the noise standard deviation and propagates it to approximate parameter uncertainties using

$\mathrm{Cov}(x) \approx \sigma^2 (A^T A)^{-1}$.

Notes:

Solving via normal equations is not very stable. 

The script assumes the third parameter corresponds to $\gamma^2$ and reports $\gamma = \sqrt{x_2}$ with propagated uncertainty.

Dependencies: `numpy`

Run: `python noisy_matrix.py`

jacobi_method_error_propagation.py

Jacobi eigenvalue algorithm for symmetric matrices. Applies plane rotations to reduce off-diagonal elements and iteratively diagonalize the matrix. After a fixed number of sweeps, the diagonal entries approximate the eigenvalues.

The script also explores eigenvalue sensitivity by perturbing a random matrix entry and measuring

$\dfrac{|\Delta \lambda|_2}{|\lambda|_2,|\delta|}$,

where $\lambda$ denotes the eigenvalue vector and $\delta$ is the relative perturbation magnitude applied to one matrix element.

Dependencies: `numpy`, `scipy`.

Run: `python jacobi_method_error_propagation.py`

### `companion_matrix_roots.f90`

Computes the roots of a polynomial by converting it into an eigenvalue problem.

Given coefficients (ascending order)

$P(x) = p_0 + p_1 x + \dots + p_n x^n$,

construct a companion matrix $C$ such that its eigenvalues are the polynomial roots:

${x: P(x)=0} = {\lambda_i(C)}$.

The program calls LAPACK DGEEV to compute eigenvalues and reports them as complex roots.

Notes:

Polynomial root finding can be ill-conditioned: small coefficient perturbations can move roots significantly (especially for high degree or clustered roots).

Requires LAPACK/BLAS.

Build/run:
`gfortran companion_matrix_roots.f90 -O2 -llapack -lblas -o companion_matrix_roots && ./companion_matrix_roots`
