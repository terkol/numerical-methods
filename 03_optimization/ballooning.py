import numpy as np
import matplotlib as mpl 
import matplotlib.pyplot as plt
from scipy.sparse import bmat  
from scipy.sparse.linalg import eigs 

def set_plot_defaults():  # Make the plots look nice
    mpl.rcParams.update({"axes.grid": True, 
                        "axes.linewidth": 0.2,  
                        "grid.linewidth": 0.2,
                        "lines.linewidth": 1, 
                        "font.family": "serif",
                        "font.serif": ["Times New Roman"], 
                        "font.size": 10,
                        "text.color": "dimgray", 
                        "xtick.color": "dimgray",
                        "ytick.color": "dimgray", 
                        "axes.labelcolor": "dimgray"
                        })

def solve_analytical(ky, g, kx=1): 
    """
    Compute the two roots of the uncoupled dispersion relation analytically.
    """
    k2 = kx**2+ky**2  
    a = 1
    b = 1j + 1j/k2 
    c = -1/k2 + g*ky*(1j + ky)/k2 
    return np.roots([a, b, c]) 

def plot_analytical(g_arr, ky_min=0, ky_max=2, N=200):  
    """
    Overlay real and imaginary parts of analytical roots for different g values.
    """
    ky_arr = np.linspace(ky_min, ky_max, N) 
    for g in g_arr:
        w_all = np.array([solve_analytical(ky, g) for ky in ky_arr]) 
        w0, w1 = w_all[:, 0], w_all[:, 1] 
        plt.plot(ky_arr, w0.real, label=f"Re ω, g={g}, root 1 (ana)")   
        plt.plot(ky_arr, w0.imag, "--", label=f"Im ω, g={g}, root 1 (ana)")
        plt.plot(ky_arr, w1.real, label=f"Re ω, g={g}, root 2 (ana)")  
        plt.plot(ky_arr, w1.imag, "--", label=f"Im ω, g={g}, root 2 (ana)")
    plt.xlabel(r"$k_y$")
    plt.ylabel(r"$\omega$")
    plt.title("Analytical eigenvalues")
    plt.legend(ncol=3)

def solve_numerical(g, M, r0, kx=1, Ln=1):
    """
    Solve the uncoupled 2x2 problem numerically for a series of ky values.
    """
    ky_arr = np.arange(1, M+1)/r0 
    mats = [np.array([[-1j, 1j+ky/Ln], [-g*ky/(kx**2+ky**2), -1j/(kx**2+ky**2)]])
        for ky in ky_arr] 
    w = np.zeros((M, 2), dtype=complex) 
    for i, A in enumerate(mats): 
        w[i, :] = np.linalg.eigvals(A) 
    return w

def plot_numerical(g_arr, M, r0, ky_min=0, ky_max=2, N=128, conv=False):
    """
    Overlay numerical roots on the analytical plot, optionally debug convergence.
    """
    ky_arr = np.linspace(ky_min, N/r0, N)      # ky vector for plotting
    for g in g_arr:
        w_all = solve_numerical(g, M, r0)[:N]  # get first N modes
        w0, w1 = w_all[:, 0], w_all[:, 1]    
        plt.plot(ky_arr, w0.real, label=f"Re ω, g={g}, root1 (num)") 
        plt.plot(ky_arr, w0.imag, "--", label=f"Im ω, g={g}, root1 (num)")
        plt.plot(ky_arr, w1.real, label=f"Re ω, g={g}, root2 (num)")  
        plt.plot(ky_arr, w1.imag, "--", label=f"Im ω, g={g}, root2 (num)")
        if g == 0.5 and conv:
            print("Convergence check, last Im ω:", w1[127].imag) 
    plt.xlabel(r"$k_y$") 
    plt.ylabel(r"$\omega$")  
    plt.title("Numerical eigenvalues vs analytical ones") 
    plt.legend(ncol=3)

def solve_coupled(g, M, r0, kx=1, Ln=1):
    """
    Build and solve the full 2M x 2M tridiagonal-block operator.
    """
    ky = np.arange(1, M+1)/r0
    layout = [[None]*M for _ in range(M)]
    for i, k in enumerate(ky):
        k2 = kx**2+k**2  
        layout[i][i] = np.array([[-1j, 1j + k/Ln], [1j - g*k/k2, -1/k2]])
        if i < M-1: 
            kp, k2p = ky[i+1], kx**2+ky[i+1]**2   # Next ky
            coeff_p = -(g*(-kx+1j*kp)) / (2*k2p)
            layout[i][i+1] = np.array([[0,0],[0,coeff_p]])  # Upper off-diagonal
        if i > 0:
            km, k2m = ky[i-1], kx**2+ky[i-1]**2   # Previous ky
            coeff_m = -(g*(-kx+1j*km))/(2*k2m)
            layout[i][i-1] = np.array([[0,0],[0,coeff_m]])  # Lower off-diagonal
    L = bmat(layout, format="csr")  # Assemble sparse matrix
    vals, vecs = eigs(L, k=1, which="LI")  # Find eigenmode with largest Im
    return vals[0], vecs[:, 0]

def max_uncoupled_mode(g, M, r0, kx=1, Ln=1): 
    """
    Loop over uncoupled 2x2 matrices to pick the mode with max Im(ω).
    """
    best = {"im": -np.inf}
    for i in range(1, M+1):     # Scan each ky
        ky = i/r0
        A = np.array([[-1j, 1j+ky/Ln], [-g*ky/(kx**2+ky**2), -1j/(kx**2+ky**2)]])
        w, v = np.linalg.eig(A) # Diagonalize
        for j in range(2):      # Check both roots
            if w[j].imag > best["im"]:  # Compare growth rate
                best.update({"im": w[j].imag, "w": w[j], "v": v[:, j], "ky": ky, "idx": i-1})  # Store best imaginary component
    return best 

def plot_mode_density(psi, ky, tag):
    """
    Display the poloidal mode envelope and its real-space transform.
    """
    M = len(ky)  
    nhat = psi[::2]                 # Extract density components
    full = np.zeros(2*M, dtype=complex)  # Prepare full spectrum
    full[:M] = nhat                 # Fill first half
    full[M:] = np.conj(nhat[::-1])  # Mirrored part for realness
    y = np.linspace(0, 2*np.pi, 2*M, endpoint=False)  # Real-space grid
    n_y = np.fft.ifft(full) * (2*M) # Inverse FFT back to y space
    plt.figure(figsize=(10, 4))
    plt.subplot(1, 2, 1) 
    plt.plot(ky, np.abs(nhat))
    plt.xlabel(r"$k_y$")
    plt.title(f"|n(k_y)| — {tag}")  
    plt.subplot(1, 2, 2)
    plt.plot(y, n_y.real) 
    plt.xlabel("y")
    plt.title(f"n(y) — {tag}")

if __name__ == "__main__": 
    set_plot_defaults()  
    M = 1024 
    r0 = 64   
    plt.figure(figsize=(12,6))  
    plot_analytical([0.0, 0.5, -0.5]) 

    plt.figure(figsize=(12,6)) 
    plot_analytical([0.0, 0.5, -0.5])
    plot_numerical([0.0, 0.5, -0.5], M, r0) 

    w_c, psi_c = solve_coupled(0.5, M, r0)  
    best = max_uncoupled_mode(0.5, M, r0)  
    best["v"] /= np.linalg.norm(best["v"])  # Normalize vector
    psi_u = np.zeros(2*M, dtype=complex)    # Initialize uncoupled psi
    psi_u[best["idx"]*2:best["idx"]*2+2] = best["v"] # Insert best mode
    ky_arr = np.arange(1, M+1) / r0         # Recompute ky grid
    plot_mode_density(psi_u, ky_arr, "uncoupled")
    plot_mode_density(psi_c, ky_arr, "coupled") 
    plt.show()  