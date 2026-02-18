import numpy as np
from scipy.optimize import minimize

def f(vars):
    x, y = vars
    return 10*np.exp(-(x-3)**2)+x**2+x*y+y**2-5*np.exp(-(y+2)**2)

if __name__ == "__main__":
    print("Nelder-Mead (0,0):")
    init = np.array([0, 0])
    nelder_mead = minimize(f, init, method='Nelder-Mead', options={'disp':True})
    print("Nelder-Mead (5,-5):")
    init = np.array([5, -5])
    nelder_mead = minimize(f, init, method='Nelder-Mead', options={'disp':True})

    print("Powell (0,0):")
    init = np.array([0,0])
    powell = minimize(f, init, method='Powell', options={'disp':True})
    print("Powell (5,-5):")
    init = np.array([5,-5])
    powell = minimize(f, init, method='Powell', options={'disp':True})
