import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

def f(x): # Given function
    return 1 / (1 + x**2)

def lagrange_method(points, y):
    pol = 0 
    for i in range(len(points)):
        lag = 1
        for j in range(len(points)):
            if i != j:
                lag = lag*((x-points[j])/(points[i]-points[j])) # Calculate langrangian from loop indexes but only when i != j
        pol += y[i] * lag # Assemble polynomial
    return sp.simplify(pol) # Return simplified

def newton_method(p, y):
    pol = y[0] 
    pol += ((y[1]-y[0])/(p[1]-p[0]))*(x-p[0]) # First divided difference
    pol += (((y[2]-y[1])/(p[2]-p[1]))-((y[1]-y[0])/(p[1]-p[0])))/(p[2]-p[1])*(x-p[0])*(x-p[1]) # Second divided difference, good thing a third wasn't required
    return sp.simplify(pol)

if __name__ == "__main__":
    points = [0, 1, 3] # Array of given points
    y = [f(x) for x in points]  # Values of y for points
    x = sp.Symbol('x') # Declare x as a symbol

    lagrange = lagrange_method(points, y)
    newton = newton_method(points, y)
    interval = np.linspace(-0.5, 3.5)
    lagrange_func = sp.lambdify(x, lagrange) # Make polynomial into function
    fx = [f(x) for x in interval]   # Values of f(x) for plot
    pol = [lagrange_func(x) for x in interval]  # Values of polynomial function

    if sp.simplify(lagrange - newton).subs(x,1) < 1e-10: # Check if polynomials are close to each other
        print("Polynomials are basically equal")

    plt.plot(interval, fx, label="Function f(x)")
    plt.plot(interval, pol, label="Interpolated polynomial")
    plt.scatter(points, y, label="Points")
    plt.legend()
    plt.grid()
    plt.show()