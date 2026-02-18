import numpy as np
import matplotlib.pyplot as plt

def f(y, z):
    return np.array([z, -y])    # returns [y', z'] = [z, -y]

def runge_kutta_solver(a, b, N, y0, yp0):
    h = (b-a)/(N-1)
    ys = np.zeros(N)
    zs = np.zeros(N)
    ys[0] = y0
    zs[0] = yp0

    for i in range(N - 1): 
        yi, zi = ys[i], zs[i] 
        k1 = f(yi, zi)*h        # Slope at the beginning of the interval
        k2 = f(yi+0.5*k1[0], zi+0.5*k1[1])*h    # Slope at midpoint 
        k3 = f(yi+0.5*k2[0], zi+0.5*k2[1])*h    # Another slope at midpoint 
        k4 = f(yi+k3[0], zi+k3[1])*h            # Slope at the end
        delta = (k1+2*k2+2*k3+k4)/6            
        ys[i+1] = yi+delta[0]
        zs[i+1] = zi+delta[1]
    return ys

if __name__ == "__main__":
    x_cos = np.linspace(0, 10, 500)
    y_cos = np.cos(x_cos) # Analytical solution for reference
    plt.plot(x_cos, y_cos, label="cos(x)", linewidth=10)

    N = 300
    ys = runge_kutta_solver(0, 10, N, 1, 0)
    xs = np.linspace(0,10,len(ys))
    plt.plot(xs, ys, label=f"RK4, N={N}", c="r")

    plt.xlabel("x")
    plt.ylabel("y")
    plt.legend()
    plt.show()