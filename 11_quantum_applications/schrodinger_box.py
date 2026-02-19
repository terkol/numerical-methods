import numpy as np                                    
import matplotlib.pyplot as plt                         

def potential_in_box(N, V0, a=-1, b=1):               
    x = np.linspace(a, b, N+2)  
    h = (b-a)/(N+1)
    D = 2*np.eye(N)-np.eye(N, k=1)-np.eye(N, k=-1)  # 3-point stencil  
    V = V0*(1-x[1:-1]**2)       # Calculate potential
    D = D/(2*h*h)+np.diag(V)    # Add potential to D
    val, vec = np.linalg.eigh(D)
    M = np.zeros((N+2, N))
    M[1:-1] = vec
    return val, M, x 

def particle_in_box(N, a=-1, b=1):
    x = np.linspace(a, b, N+2)          
    h = (b-a)/(N+1)                     
    D = 2*np.eye(N)-np.eye(N, k=1)-np.eye(N, k=-1)  # Stencil  
    val, vec = np.linalg.eigh(D/(2*h*h))
    M = np.zeros((N+2, N))
    M[1:-1] = vec 
    return val, M, x                               

if __name__ == "__main__":

    # Unperturbed box  
    N = 800                                             
    val, vec, x = potential_in_box(N, V0=0)
    vec[:,0] /= np.sqrt(np.trapezoid(vec[:,0]**2, x))
    E0 = val[0]
    integral = np.trapezoid(vec[:,0]**2*(1-x**2), x)

    # Perturbations
    for V0 in [1, 10, 30]:  # Loop over strengths
        val, psi, x = potential_in_box(N, V0)
        E = val[0] 
        E_pt = E0 + V0*integral
        print(f"V0={V0:2d}: E={E:.6f}, dE={E_pt-E:.6f}")

    val, vec, x = particle_in_box(N) 

    for n in range(5):
        plt.plot(x, vec[:,n], label=f"n={n+1}, E={val[n]:.6f}")
    plt.xlabel("x")
    plt.ylabel("psi(x)")
    plt.legend()  
    plt.grid()
    plt.show()
