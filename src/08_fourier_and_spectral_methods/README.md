# Fourier and spectral methods

This folder contains small programs illustrating discrete Fourier transforms (DFT/FFT), time–frequency analysis, and spectral time stepping for a diffusion equation on a periodic domain.

## Programs
### `fourier_transform.py`

Loads a sampled time series from sample1.dat and analyzes its frequency content using an FFT and a spectrogram.

Given samples $s_n = s(t_n)$ with sampling interval $\Delta t$, the discrete Fourier transform resolves frequencies

$f_k = \dfrac{k}{N\Delta t}$ (via numpy.fft.fftfreq),

and the dominant (fundamental) frequency can be estimated by locating the largest peak in the magnitude spectrum $|\hat s(f)|$ after removing the mean (to suppress the DC component).

The script also computes a spectrogram (short-time Fourier transform) to visualize how spectral content evolves over time.

Dependencies: `numpy`, `matplotlib`, `scipy`, `pyfftw`.

Run: `python fourier_transform.py`

Data:

expects `sample1.dat` in the data folder (two columns: time, signal).

### `inverse_fourier_transform.py`

Spectral solution of a diffusion/heat-type evolution on a periodic domain using FFTs.

Starting from an initial condition $u(x,0)$ sampled on a uniform grid, the script computes Fourier coefficients $\hat u_k(0)$ with fft, evolves them in time by multiplying each mode by an exponential decay factor, and reconstructs $u(x,t)$ with ifft.

For the heat equation on a periodic domain,

$u_t = u_{xx}$,

the Fourier modes evolve as

$\hat u_k(t) = \hat u_k(0),e^{-k^2 t}$,

so high-frequency modes decay faster, producing smoothing as $t$ increases.

Dependencies: `numpy`, `matplotlib`.

Run: `python inverse_fourier_transform.py`