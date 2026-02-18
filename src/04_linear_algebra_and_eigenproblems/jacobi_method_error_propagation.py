import numpy as np
import scipy.linalg as la
import random 

def jacobi(A, N, imax=10):
    A = A.copy()    
    for _ in range(imax):
        B = A.copy()
        for p in range(N-1):
            for q in range(p+1, N):
                # Make rotation matrix
                angle = 0.5 * np.arctan2(2 * B[p, q], (B[p, p] - B[q, q]))  
                c = np.cos(angle)
                s = np.sin(angle)
                Q = np.eye(N)
                Q[p, p] = c
                Q[q, q] = c
                Q[p, q] = -s
                Q[q, p] = s
                B = Q.T @ B @ Q # Zero out off-diagonal element using Q
        A = B
    return np.diag(A)

if __name__ == "__main__":
    Q = np.array([[4.2, 4.6, 3.9],  [4.6, 3.3, 2.7], [3.9, 2.7, 2.1]]) # Test matrix
    print("Jacobi: ", jacobi(Q, len(Q)))
    print("Scipy:  ", la.eigvalsh(Q))

    def err_propag(N, delta):
        Q = np.eye(N)
        for p in range(N):
            for q in range(N): # Reassigning half of values, I know, but it looks neater 
                rand = random.random()
                Q[p, q] = rand  
                Q[q, p] = rand  
        eig = jacobi(Q, N)  # Eigenvalues
        eig_norm = sum(eig**2)**0.5 # Norm of eigenvalues
        randp = int(random.random()*N)
        randq = int(random.random()*N)
        Q[randp, randq] = Q[randp, randq] * (1.0 + delta)
        Q[randq, randp] = Q[randp, randq]
        eigp = jacobi(Q, N) # Eigenvalues of perturbed matrix
        err = eigp - eig
        err_norm = sum(err**2)**0.5 # Norm of the difference
        return err_norm/(eig_norm*delta)

    f = []
    # Nested loop cycling through different system sizes and perturbations
    for N in range(8, 12):
        for q in range(20):
            delta = 10.0**(-(q+1))     # 1e-1, 1e-2, ..., 1e-20 (pick range you like)
            f.append(err_propag(N, delta))
    
    f = np.array(f)
    fmean = sum(f)/len(f)
    print("Mean:", fmean)
    print("SD:  ", (sum((f-(fmean))**2)/len(f))**0.5)