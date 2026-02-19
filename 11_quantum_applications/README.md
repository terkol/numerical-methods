# Quantum applications

This folder contains small programs that apply numerical linear algebra and finite-difference discretization to quantum-mechanical eigenvalue problems.

## Programs

### `schrodingers_box.py`

Finite-difference solution of the 1D time-independent Schrödinger equation on a bounded interval with Dirichlet boundary conditions.

In units where $\hbar = 1$ and $m = 1$, the Hamiltonian is

$H = -\dfrac{1}{2}\dfrac{d^2}{dx^2} + V(x)$,

and stationary states satisfy the eigenvalue problem

$H\psi = E\psi$,

with $\psi(a)=\psi(b)=0$ (“particle in a box” boundary conditions).

Discretization:

The second derivative is approximated on the interior grid points using the standard 3-point stencil:

$-\dfrac{d^2\psi}{dx^2}(x_i) \approx \dfrac{2\psi_i - \psi_{i+1} - \psi_{i-1}}{h^2}$,

so the discrete operator is

$H_h = \dfrac{1}{2h^2}(2I - T_{+1} - T_{-1}) + \mathrm{diag}(V)$.

The script solves two cases:

Unperturbed box ($V(x)=0$), computing low-lying eigenvalues/eigenfunctions.

Quadratic perturbation $V(x)=V_0(1-x^2)$, comparing the numerical ground-state energy shift to first-order perturbation theory:

$E_0(V_0) \approx E_0(0) + V_0 \int |\psi_0(x)|^2(1-x^2),dx.$

Outputs:

Prints numerical ground-state energies for several $V_0$ and the difference from the first-order perturbation estimate.

Plots the first few eigenfunctions for the unperturbed box.

Dependencies: `numpy`, `matplotlib`.