import numpy as np

def f(y, z):
    return np.array([z, -y])

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

def shooting_solver(a, b, N, y0, y1, ka=0, kb=2, max_i=30, tol=1e-6):
    def g(k):
        return runge_kutta_solver(a, b, N, y0, k)[-1] - y1  # Boundary error at x=b  
    if g(ka)*g(kb) < 0:
        print('No bracket')
    g_a = g(ka)                 # Error at low k  
    for _ in range(max_i):      # Bisection iterations  
        km = 0.5*(ka + kb)      # Midpoint k  
        g_m = g(km)             # Error at midpoint  
        if g_a * g_m < 0:       # Sign change on [ka, km]?  
            kb = km             # Bracket right side  
        else:
            ka, g_a = km, g_m   # Bracket left side  
        if abs(g_m) < tol:
            break
    k_star = 0.5*(ka + kb)  # Final root estimate  
    return k_star               

if __name__ == "__main__":
    for i in [10, 30, 100]:
        print(f"N={i}, k* = {shooting_solver(0, 1, i, 1, 1)}")
