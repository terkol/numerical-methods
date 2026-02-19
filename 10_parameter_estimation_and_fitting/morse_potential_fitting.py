import numpy as np
import matplotlib.pyplot as plt
import os

def f(r, x):
    return x[0]*(1-np.exp(-x[1]*(r-x[2])))**2   # Morse potential for given r and x

def obj(r, V, x):
    return sum([0.5*(f(r[i], x)-V[i])**2 for i in range(len(r))])   # Objective function for multiple r and V

def ddx1(r, x):
    return (1-np.exp(-x[1]*(r-x[2])))**2 # Partial derivative in regards to x1

def ddx2(r, x):
    e = np.exp(-x[1]*(r-x[2]))
    return 2*x[0]*(1 - e)*(r - x[2])*e    # Partial derivative in regards to x2

def ddx3(r, x):
    return -2*x[0]*x[1]*(1-np.exp(-x[1]*(r-x[2])))*np.exp(-x[1]*(r-x[2]))       # Partial derivative in regards to x3

def ddx(r, x, n):   # Select derivative with n
    if n == 0:
        return ddx1(r, x)
    elif n == 1:
        return ddx2(r, x)
    else:
        return ddx3(r, x)

def jacobian(r, x): # Assemble Jacobian
    J = []
    for i in range(len(r)): # For every r in r_data...
        J.append([ddx(r[i], x, j) for j in range(len(x))])  # ...Append a row of derivatives (x1, x2, x3) into J with that r
    return np.array(J)

def f_v(r, x):
    return np.array([x[0]*(1-np.exp(-x[1]*(i-x[2])))**2 for i in r])  # Construct f(x) vector with r_data

def gradient(x, y): # Calculate gradient
    r = r_data
    J = jacobian(r, x)
    fxy = f_v(r, x)-y
    return J.T @ fxy

def gradient_descent(r, x, gamma, y, max_i=22000):  # Move into the direction of negative descent by gamma amount
    objs = []
    sum = 1
    for i in range(max_i):
        x = x - gamma*gradient(x, y)
        if i%100 == 0:              # Every 100 steps, count new value for objective function and add it into an array
            sum = obj(r, y, x)
            objs.append(sum)
    return x, objs                  # Return final x and list of objective function values
    
if __name__ == "__main__":
    data = np.loadtxt(os.path.dirname(os.path.dirname(__file__))+'\\data\\raw\\morse_data.txt')
    r_data = data[:,0]
    V_data = data[:,1]
    plt.scatter(r_data, V_data)
    plt.grid()
    plt.show()
    print(f"Sum of squares: {obj(r_data, V_data, (0.33, 2, 2.8))}")
    x, objs = gradient_descent(r_data, (0.33, 2, 2.8), 1e-4, V_data)

    print(f"Final x: {x}")
    print(f"final objective function value: {objs[-1]}")

    plt.scatter(r_data, V_data, label="Data")
    plt.plot(np.linspace(1.7, 10, 200),np.array([f(ri, x) for ri in np.linspace(1.7, 10, 200)]), label="Fitted Potential")
    plt.xlabel("r")
    plt.ylabel("V")
    plt.legend()
    plt.grid()
    plt.show()

    plt.plot(list(range(0,22000,100)),objs)
    plt.xlabel("Iteration")
    plt.ylabel("Objective Function")
    plt.grid()
    plt.show()