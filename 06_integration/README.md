# Numerical integration

This folder contains small programs for numerical quadrature. Topics include Romberg integration (Richardson extrapolation of the trapezoidal rule), composite Simpson refinement, and Gauss–Laguerre quadrature for integrals on $[0,\infty)$.

## Programs

### `romberg_method.py`

Romberg integration for approximating a definite integral

$I = \int_a^b f(x),dx$.

Romberg builds a table $R_{n,m}$ using trapezoidal-rule refinements and Richardson extrapolation. The first column is the composite trapezoid rule with step size $h/2^n$:

$R_{0,0} = \dfrac{h}{2},[f(a)+f(b)]$, where $h=b-a$,

and successive refinements add evaluations at the new midpoints. Higher columns apply Richardson extrapolation assuming an even-power error expansion:

$R_{n,m} = \dfrac{4^m R_{n,m-1} - R_{n-1,m-1}}{4^m - 1}$.

The script tests two integrands on $[0,1]$ involving $\sin(\sqrt{x})$. In the second case it integrates

$\sin(\sqrt{x}) - \sqrt{x}$

and adds back the known analytic integral

$\int_0^1 \sqrt{x},dx = 2/3$,

which can improve numerical behavior near $x=0$.

Dependencies: `numpy`.

Run: `python romberg_method.py`

### `composite_simpson_refinement.py`

Composite Simpson refinement and Gauss–Laguerre quadrature examples.

Simpson’s rule on an interval $[a,b]$ is

$S(a,b) = \dfrac{b-a}{6},\Big(f(a) + 4f\big(\tfrac{a+b}{2}\big) + f(b)\Big)$.

The script repeatedly refines the interval by splitting into $2^k$ equal subintervals and summing Simpson’s rule over each subinterval until successive refinements change by less than a tolerance (a simple convergence heuristic).

It also demonstrates Gauss–Laguerre quadrature for integrals of the form

$\int_0^\infty e^{-x} g(x),dx \approx \sum_{i=1}^n w_i, g(x_i)$,

using `numpy.polynomial.laguerre.laggauss(n)` to obtain nodes $x_i$ and weights $w_i$.

Notes:

The Simpson refinement in this script is uniform (not locally adaptive).

The change between successive refinements is a heuristic error indicator; for Simpson’s method a common scaling is $|S_{2h}-S_h|/15$.

The variable substitution $x = t/(1-t)$ maps $[0,1)$ to $[0,\infty)$ and introduces a singular limit as $t\to 1$; the script truncates at $t=0.999$.

Dependencies: `numpy`.

Run: `python composite_simpson_refinement.py`