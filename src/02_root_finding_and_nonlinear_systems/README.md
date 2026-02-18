# Root finding and nonlinear systems

This folder contains small programs for solving nonlinear equations and nonlinear systems. It includes bracketing methods (bisection), locally convergent methods (Newton), quasi-Newton updates (Broyden), and a polynomial-root example computed via an eigenvalue (companion matrix) formulation.

## Programs

### `bisection_vs_newton.f90`

Compares bisection and Newton’s method on a scalar equation $f(x)=0$.

Bisection maintains an interval $[a,b]$ such that

$f(a),f(b) < 0$,

and repeatedly halves it by evaluating $f(c)$ at

$c = (a+b)/2$,

keeping the subinterval that preserves the sign change. For continuous $f$, this converges globally to a root in the bracket.

Newton’s method updates

$x_{k+1} = x_k - \frac{f(x_k)}{f'(x_k)}$,

and converges quickly when started sufficiently close to a simple root (typically quadratic convergence). In this implementation, $f'(x)$ is approximated by a central difference:

$f'(x) \approx \dfrac{f(x+h)-f(x-h)}{2h}$.

Notes:

Bisection requires $f$ to be continuous on $[a,b]$ and a valid sign-change bracket.

The specific $f(x)$ used in the file has singularities at $x=\pm 1$ due to a division by $(x^2-1)$, so brackets/iterates should avoid those points.

Build/run:
`gfortran bisection_vs_newton.f90 -O2 -o bisection_vs_newton && ./bisection_vs_newton`

### `broydens_method.py`

Broyden’s method (a quasi-Newton method) for a nonlinear system in $\mathbb{R}^2$:

$F(x) = 0$, where $x = (x_1,x_2)$.

The iteration uses an approximate inverse Jacobian $B_k \approx J(x_k)^{-1}$:

$x_{k+1} = x_k - B_k,F(x_k)$.

After each step, a rank-1 update enforces a secant condition using

$s_k = x_{k+1}-x_k$,
$y_k = F(x_{k+1})-F(x_k)$,

and updates $B_k$ instead of recomputing a Jacobian.

The example system is:

$f_1(x) = \exp(-(x_1^2 + x_2^2)) - \frac18$
$f_2(x) = \sin(x_1) - \cos(x_2)$.

Notes:

Convergence depends strongly on the initial guess.

The rank-1 update can become unstable if the update denominator becomes very small (too lazy to add safeguards).

Dependencies: `numpy`

Run: `python broydens_method.py`

### `find_root_eigenvalue.py`

Polynomial roots computed as eigenvalues of a companion matrix.

Given coefficients in ascending order

$P(x) = p_0 + p_1 x + \dots + p_n x^n$ (with $p_n \neq 0$),

one can form the companion matrix $C \in \mathbb{R}^{n\times n}$ such that its eigenvalues are the roots of $P$:

${\lambda_i(C)} = {x : P(x)=0}$.

This script compares the eigenvalue-based roots to numpy.roots (which expects coefficients in descending order).

Dependencies: `numpy`, `scipy`.

Run: `python find_root_eigenvalue.py`

### `logistic_fixed_point_function_plot.py`

Plots the function

$f(x;\mu) = x - \mu x(1-x)$

for several values of $\mu$. This function arises when studying fixed points of the logistic map $x \mapsto \mu x(1-x)$ by rewriting the fixed-point condition as a root-finding problem:

$x = \mu x(1-x)\quad \Longleftrightarrow \quad f(x;\mu)=0$.

Dependencies: `numpy`, `matplotlib`.

Run: `python logistic_fixed_point_function_plot.py`