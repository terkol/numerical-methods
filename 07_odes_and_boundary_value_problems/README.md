# ODEs and boundary value problems

This folder contains small programs for ordinary differential equations (ODEs), including initial value problems (IVPs) solved with 4th-order Runge–Kutta (RK4) and a boundary value problem (BVP) solved by the shooting method.

## Programs
### `rk4.py`

Fourth-order Runge–Kutta solver for the second-order ODE

$y'' + y = 0$,

rewritten as a first-order system with $z = y'$:

$y' = z$

$z' = -y$.

Given initial conditions $(y(a), y'(a)) = (y_0, z_0)$, RK4 advances the state with step size $h$ using

$k_1 = h,f(y_n, z_n)$
$k_2 = h,f(y_n + k_1/2, z_n + k_1/2)$
$k_3 = h,f(y_n + k_2/2, z_n + k_2/2)$
$k_4 = h,f(y_n + k_3, z_n + k_3)$

and updates with the weighted average

$(y_{n+1}, z_{n+1}) = (y_n, z_n) + (k_1 + 2k_2 + 2k_3 + k_4)/6$.

The script compares the numerical solution to the analytic solution $y(x)=\cos(x)$ for $(y(0),y'(0))=(1,0)$.

Dependencies: `numpy`, `matplotlib`.

Run: `python rk4.py`

### `shooter.py`

Shooting method for a two-point boundary value problem (BVP) based on the same ODE

$y'' + y = 0$,

with boundary conditions

$y(a) = y_0$, 

$y(b) = y_1$.

The method converts the BVP into a root-finding problem for the unknown initial slope $k = y'(a)$:

Define the boundary mismatch

$g(k) = y(b; k) - y_1$,

where $y(b;k)$ is obtained by solving the IVP with RK4 starting from $(y(a),y'(a))=(y_0,k)$. The correct slope satisfies

$g(k^*) = 0$.

This implementation uses bisection on an initial bracket $[k_a, k_b]$ to find $k^*$.

Dependencies: `numpy`.

Run: `python shooter.py`