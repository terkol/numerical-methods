import numpy as np  

def romberg_method(f, a, b, N=5):
    RE = np.zeros(N)
    h = b-a
    RE[0, 0] = 0.5*h*(f(a)+f(b))        # first trapezoidal approximation 
    for n in range(1, N): 
        h /= 2              # Halve interval length each loop
        s = 0                   
        for k in range(1, 2**n, 2):
            s += f(a+k*h)             # Add function value to sum
        RE[n, 0] = 0.5*RE[n-1, 0]+s*h # Update the table with the new approximation
        for m in range(1, n+1):                     
            RE[n, m] = (4**m*RE[n][m-1]-RE[n-1][m-1])/(4**m-1)  # Use Richardson's method to refine approximation
    return RE[N-1, N-1]                 # Return the most accurate approximation from the Romberg table

def f_b(x):
    return np.sin(np.sqrt(x)) 

def f_c(x):  
    return np.sin(np.sqrt(x)) - np.sqrt(x) 

if __name__ == "__main__":
    print("Integral of f_b(x):", romberg_method(f_b, 0, 1))  
    print("Integral of f_c(x):", romberg_method(f_c, 0, 1) + 2/3)  
