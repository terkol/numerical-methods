import numpy as np
from scipy.linalg import eigvals

def myroots(N, p):
    n = N - 1
    A = np.zeros((n, n))
    A[1:,:-1] = np.eye(n-1) # subdiagonal ones
    A[:,-1] = -p[:-1] / p[-1]
    return eigvals(A)

if __name__ == "__main__":
    p = [-6, 11, -6, 1] # Third degree polynomial
    N = len(p)

    roots_eigenvalues = myroots(N, p)
    print("myroots: ", np.sort(roots_eigenvalues))

    roots_np = np.roots(p[::-1])  # Need to reverse p because np.roots wants highest degree first
    print("np.roots:", np.sort(roots_np))