import numpy as np

def derivative(x=1.0):
    return (-2*x*np.exp(-x**2))*np.cos(np.exp(-x**2)) # Precise derivative at x = 1 for comparison

def func(x):
    return np.sin(np.exp(-x**2))  

def central_difference(h, x=1): 
    return (func(x+h)-func(x-h))/(2*h)

def richardson(N=5, h=0.1):
    RE = np.eye(N)  
    for n in range(N):
        RE[n][0] = central_difference(h/2**n)
    for n in range(1,N):
        for m in range(1,n+1):
            RE[n][m] = (4**m*RE[n][m-1]-RE[n-1][m-1])/(4**m-1)
    return RE[N-1, N-1] # return most precise approximation

if __name__ == "__main__":
    print(f"Analytical derivative:   {derivative()}")
    print(f"Richardson approximation: {richardson()}")

    # Print approximations with different h values using central difference
    for i in range(2,8):
        print(f"Central difference for h=10**-{i*2}: {central_difference(10**-(i*2))}")