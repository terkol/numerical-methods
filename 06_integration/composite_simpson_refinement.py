import numpy as np

def composite_simpsons(func, a, b, tol=1e-6): 
    integral = simpsons(func, a, b)
    error = 1
    subdivisions = 2 
    while abs(error) > tol: 
        intervals = np.linspace(a, b, 2**subdivisions)
        finer_integral = sum(simpsons(func, intervals[i], intervals[i+1]) for i in range(len(intervals)-1))
        error = finer_integral - integral
        integral = finer_integral
        subdivisions += 1
    return integral

def simpsons(func, a, b):
    return ((b-a)/6)*(func(a)+4*func((a+b)/2)+func(b))  

def f(x): 
    return np.exp(-x)*np.cos(x)**2 

def f_t(t): 
    return (np.exp(-(t/(1-t)))*np.cos(t/(1-t))**2)/(1-t)**2

def gauss_laguerre(g, n): 
    nodes, weights = np.polynomial.laguerre.laggauss(n)
    return np.sum(weights * g(nodes))

def g(x):
    return np.cos(x)**2  

if __name__ == "__main__":
    print("f(x):   " + str(composite_simpsons(f, 0, 1000)))  
    print("f_t(x): " + str(composite_simpsons(f_t, 0, 0.999))) 
    print("Gauss-Laguerre, n=2:  " + str(gauss_laguerre(g, 2)))
    print("Gauss-Laguerre, n=4:  " + str(gauss_laguerre(g, 4)))
    print("Gauss-Laguerre, n=8:  " + str(gauss_laguerre(g, 8)))
