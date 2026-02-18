# Parameter estimation and fitting

This folder contains small programs for estimating model parameters from data. It includes a basic Metropolis MCMC sampler for linear regression with Gaussian errors, and gradient-based least-squares fitting of a Morse potential.

## Programs
### `likelihood_fitting.py`

Metropolis (random-walk) MCMC for a straight-line model with heteroscedastic Gaussian errors:

$y_i = x_1 t_i + x_2 + \varepsilon_i,\qquad \varepsilon_i \sim \mathcal{N}(0,\sigma_i^2).$

Gaussian likelihood:

$L(x_1,x_2) = \prod_i \frac{1}{\sqrt{2\pi}\sigma_i}\exp\!\left(-\frac{(y_i-(x_1 t_i + x_2))^2}{2\sigma_i^2}\right).$

Metropolis acceptance probability:

$\alpha = \min\!\left(1,\; \frac{L(x_1',x_2')}{L(x_1,x_2)}\right).$

Numerical note (intentional instability):
- This script intentionally computes $L$ as a *direct product* of many small floating-point numbers to demonstrate numerical underflow/overflow.
- For moderately sized datasets, the product often underflows to $0$ (or occasionally overflows), making the acceptance ratio $q = L_{\text{new}}/L_{\text{old}}$ become $0$, $\infty$, or NaN (e.g. the indeterminate case $0/0$). This can distort or stall the Markov chain.
- The implementation prints simple counters (e.g. how often $L=0$ or $q$ becomes NaN) to make the failure mode visible.

Stable alternative (not used here):
- The standard fix is to work with the log-likelihood $\log L$ and accept/reject using differences $\log L_{\text{new}}-\log L_{\text{old}}$, which avoids underflow. This repository keeps the raw-product form on purpose for demonstration.

Dependencies: `numpy`, `matplotlib`.

Run: `python likelihood_fitting.py`

Data:
- expects `linear_data_errors.txt` in the same folder (columns: $t$, $y$, $\sigma$).

### `morse_potential_fitting.py`

Least-squares fitting of the Morse potential

$V(r) = D,(1-e^{-a(r-r_0)})^2$

to data $(r_i, V_i)$ by minimizing the objective

$J(D,a,r_0) = \frac12 \sum_i (V(r_i;D,a,r_0) - V_i)^2.$

The script computes the Jacobian (partial derivatives of $V$ with respect to parameters) and performs gradient descent:

$x_{k+1} = x_k - \gamma \nabla J(x_k)$,

with a fixed step size $\gamma$, while tracking the objective value every 100 iterations. It plots the fitted curve against the data and the objective value versus iteration.

Dependencies: `numpy`, `matplotlib`.

Run: `python morse_potential_fitting.py`

Data:
- expects `morse_data.txt` in the data folder (columns: $r$, $V$).