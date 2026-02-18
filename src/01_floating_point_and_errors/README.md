# Floating point and numerical error

This folder contains small programs illustrating common floating-point issues in scientific computing: loss of significance (catastrophic cancellation), rounding-error accumulation in summation, and truncation vs rounding tradeoffs in numerical differentiation.

## Programs

### `exp.py`

Stable evaluation of the exponential by range reduction. Write

$x = k \ln(2) + r$, with $k$ an integer,

so that $r$ is “small”, then compute

$exp(x) = 2^k exp(r)$.

This reduces numerical risk (overflow/underflow) and can improve accuracy when $|x|$ is large.

Run: `python exp.py`

### `harmonic_sum_naive.f90`

Naive harmonic series summation

$H_N$ = sum from $n=1$ to $N$ of $1/n$,

and a “batched” summation that changes the association of additions. Demonstrates that floating-point addition is not associative:

$(a + b) + c$ is not necessarily equal to $a + (b + c)$

in finite precision, so different summation groupings can produce different rounding error.

Build/run:
`gfortran harmonic_sum_naive.f90 -O2 -o harmonic_sum_naive ./harmonic_sum_naive`

### `harmonic_sum_kahan.f90`

Kahan compensated summation applied to the harmonic series. Kahan’s method maintains a compensation term c that tracks low-order bits lost to rounding. In one common form:

$y = x - c$

$t = s + y$

$c = (t - s) - y$

$s = t$

This typically reduces error relative to naive summation, especially in single precision.

Build/run:
`gfortran harmonic_sum_kahan.f90 -O2 -o harmonic_sum_kahan ./harmonic_sum_kahan`

### `richardson_vs_analytical.py`

Central difference approximation for the derivative:

$f'(x) \approx ( f(x+h) - f(x-h) ) / (2h)$

and Richardson extrapolation using the error model $D(h) = f'(x) + α h^2 + O(h^4)$. The extrapolation table is built by:

$R(n,0) = D(h / 2^n)$

$R(n,m) = ( 4^m R(n,m-1) - R(n-1,m-1) ) / (4^m - 1)$

The test function is:

$f(x) = \sin( \exp(-x^2) )$

with derivative:

$f'(x) = \cos( \exp(-x^2) ) * ( -2x \exp(-x^2) )$

At x = 1:

$f'(1) = -2 \exp(-1) \cos( \exp(-1) )$

This example shows the tradeoff between truncation error (large $h$) and rounding/cancellation error (very small $h$).

Run: `python richardson_vs_analytical.py`

### `taylor_stable_evolution.f90`

Stable evaluation near $x = 0$ using Taylor expansions for expressions with removable singularities / cancellation:

$f_1(x) = (\cos(x) - 1) / x^2$
$f_2(x) = (\exp(x) - \exp(-x)) / (2x) = \sinh(x) / x$

Series expansions around $x = 0$:

$\cos(x) = 1 - x^2/2 + x^4/24 + O(x^6)$
so $f_2(x) = -1/2 + x^2/24 + O(x^4)$

$\sinh(x) = x + x^3/6 + O(x^5)$
so $f_2(x) = 1 + x^2/6 + O(x^4)$

Demonstrates why series approximations are preferred for small $|x|$ to avoid catastrophic cancellation.

Build/run:
`gfortran taylor_stable_evolution.f90 -O2 -o taylor_stable ./taylor_stable`
