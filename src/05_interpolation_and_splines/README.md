# Interpolation and splines

This folder contains small programs illustrating polynomial interpolation and spline construction. Topics include Lagrange vs Newton forms of the interpolating polynomial, B-spline basis functions, and cubic spline interpolation with different boundary conditions.

## Programs

### `interpolation_methods.py`

Interpolates the function

$f(x) = \dfrac{1}{1+x^2}$

through a small set of nodes using two equivalent polynomial forms:

- Lagrange form (for nodes $x_i$ with values $y_i$):
$P(x) = \sum_{i=0}^{n-1} y_i , \ell_i(x)$,
$\ell_i(x) = \prod_{j\ne i} \dfrac{x-x_j}{x_i-x_j}$.

- Newton form (using divided differences):
$P(x) = c_0 + c_1(x-x_0) + c_2(x-x_0)(x-x_1) + \dots$

In the current script the Newton polynomial is written explicitly for 3 nodes (up to the second divided difference). The script plots $f(x)$ and the interpolating polynomial on an interval containing the nodes.

Dependencies: `numpy`, `sympy`, `matplotlib`.

Run: `python interpolation_methods.py`

### `bspline_test.py`

Visualizes B-spline basis elements of increasing degree constructed from a small knot set. Uses SciPy’s BSpline.basis_element(...) to build basis functions and plots them over their support.

Dependencies: `numpy`, `scipy`, `matplotlib`.

Run: `python bspline_test.py`

### `bspline_control_points.py`

Parametric cubic B-spline curve defined by manually chosen control points. Treats the x- and y-coordinates as separate spline coefficient sequences:

$x(t) = \sum_i c_i^x B_i(t)$,
$y(t) = \sum_i c_i^y B_i(t)$,

and plots the resulting curve together with the control polygon. Demonstrates how control points shape the spline.

Dependencies: `numpy`, `scipy`, `matplotlib`.

Run: `python bspline_control_points.py`

### `scipy_spline_test.py`

Cubic spline interpolation of random data with different boundary curvature constraints. Uses SciPy CubicSpline with boundary conditions specified as endpoint second derivatives:

$S''(x_0) = \alpha$, $S''(x_n) = \alpha$,

and compares the resulting splines for different node densities and curvature values.

Dependencies: `numpy`, `scipy`, `matplotlib`.

Run: `python scipy_spline_test.py`