# Probability and statistics

This folder contains small Monte Carlo and statistical demos: bootstrap confidence intervals, the central limit theorem (CLT) and its failure cases, and a chi-squared goodness-of-fit test for random number generators.

## Programs

### `bootstrap.py`

Bootstrap resampling to estimate uncertainty in the sample mean.

Given observations $x_1,\dots,x_n$, the bootstrap draws $B$ resamples of size $n$ (sampling with replacement) and computes the mean for each resample:

$\bar{x}^{*(b)} = \dfrac{1}{n}\sum_{i=1}^n x_i^{*(b)}$.

A percentile confidence interval for the mean is then

$[\mathrm{q}{0.025}(\bar{x}^*),; \mathrm{q}{0.975}(\bar{x}^*)]$,

where $\mathrm{q}_p$ denotes the $p$-th quantile of the bootstrap distribution.

The script generates exponential samples (mean 5) and forms a bootstrap 95% interval for the mean, and can be used to empirically test coverage (fraction of intervals containing the true mean).

Dependencies: `numpy`, `matplotlib`.

Run: `python bootstrap.py`

### `central_limit_theorem.py`

Demonstrates the CLT by sampling i.i.d. random variables and plotting the distribution of sample means

$Z_K = \dfrac{1}{K}\sum_{i=1}^K X_i$,

for increasing $K$. For distributions with finite variance, the CLT predicts that (after appropriate centering and scaling) the distribution approaches a normal distribution as $K \to \infty$.

The script includes a heavy-tailed case (Cauchy/Lorentz) where the mean/variance are not finite; in that case the CLT does not apply and the distribution of sample means does not become normal.

Dependencies: `numpy`, `matplotlib`.

Run: `python central_limit_theorem.py`

`chi_squared.py`

Chi-squared goodness-of-fit test for uniformity of binned random samples.

For counts $O_i$ in $k$ bins and expected counts $E_i$ (uniform: $E_i = N/k$), Pearson’s chi-squared statistic is

$\chi^2 = \sum_{i=1}^k \dfrac{(O_i - E_i)^2}{E_i}$.

Under the null hypothesis (independent uniform draws) and for sufficiently large counts, $\chi^2$ is approximately distributed as a chi-squared random variable with

$\mathrm{df} = k-1$

degrees of freedom. The script repeatedly generates binned samples, computes $\chi^2$, and compares the histogram to the theoretical $\chi^2_{\mathrm{df}}$ density. An intentionally poor linear congruential generator (LCG) can be substituted to show deviations from the theoretical distribution.

Dependencies: `numpy`, `matplotlib`, `scipy`.

Run: `python chi_squared.py`