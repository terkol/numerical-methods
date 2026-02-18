import numpy as np

def broyden(x_init, i_max=100, tol=1e-8):
    x = np.array(x_init) 
    B = np.eye(len(x))      # Identity matrix as initial guess for inverted Jacobian
    for _ in range(i_max): 
        f = np.array([f1(x), f2(x)])    # Current correction values

        if np.linalg.norm(f) < tol: 
            break

        x_last = x.copy() 
        x = x - B @ f  
        s = x - x_last  
        y = np.array([f1(x), f2(x)]) - np.array([f1(x_last), f2(x_last)])
        denom = (s.T @ B @ y)
        if abs(denom) > 1e-10:
            B = B + np.outer((s - B @ y), s) / denom     # Update inverted Jacobian
    return x  

def f1(x):
    return np.exp(-(x[0]**2 + x[1]**2)) - 1/8

def f2(x): 
    return np.sin(x[0]) - np.cos(x[1])

if __name__ == "__main__":
    print("Broyden:", broyden([1.0, 1.0])) 