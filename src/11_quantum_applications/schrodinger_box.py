import numpy as np                                    
import matplotlib.pyplot as plt                         

def potential_in_box(N, V0, a=-1, b=1):               
    x = np.linspace(a, b, N+2)  
    h = (b-a)/(N+1)             # Stepsize
    D = 2*np.eye(N)-np.eye(N, k=1)-np.eye(N, k=-1)  # 3-point stencil  
    V = V0*(1-x[1:-1]**2)       # Calculate potential
    D = D/(2*h*h)+np.diag(V)    # Add potential to D
    val, vec = np.linalg.eigh(D)# Find eigenvalues and -vectors
    M = np.zeros((N+2, N))      # Allocate full taller matrix 
    M[1:-1] = vec               # Insert eigenvectors in the middle  
    return val, M, x            # Return eigenvals, eigenfuncs, grid  

def particle_in_box(N, a=-1, b=1):
    x = np.linspace(a, b, N+2)          
    h = (b-a)/(N+1)                     # Stepsize
    D = 2*np.eye(N)-np.eye(N, k=1)-np.eye(N, k=-1)  # Stencil  
    val, vec = np.linalg.eigh(D/(2*h*h))# Find eigenvalues and -vectors
    M = np.zeros((N+2, N))              # Allocate full taller matrix
    M[1:-1] = vec                       # Insert eigenvectors in the middle  
    return val, M, x                               

if __name__ == "__main__":
    N = 800                                             
    val, vec, x = potential_in_box(N, V0=0)             # Unperturbed box  
    vec[:,0] /= np.sqrt(np.trapezoid(vec[:,0]**2, x))   # Normalize ground state
    E0 = val[0]                                         # Unperturbed ground state energy  
    integral = np.trapezoid(vec[:,0]**2*(1-x**2), x)    # Compute unperturbed integral

    for V0 in [1, 10, 30]:                              # Loop over strengths  
        val, psi, x = potential_in_box(N, V0)           # Solve with potential  
        E = val[0]                                      # Numeric ground-state E  
        E_pt = E0 + V0*integral                         # First-order PT estimate  
        print(f"V0={V0:2d}: E={E:.6f}, dE={E_pt-E:.6f}")# Report numeric vs PT  

    val, vec, x = particle_in_box(N) 

    for n in range(5):
        plt.plot(x, vec[:,n], label=f"n={n+1}, E={val[n]:.6f}")
    plt.xlabel("x")
    plt.ylabel("psi(x)")
    plt.legend()  
    plt.grid()
    plt.show()
