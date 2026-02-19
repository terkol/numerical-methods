import numpy as np
import matplotlib.pyplot as plt

def f(x):
    return (x**2)*(pos_erf(x)-1)

def pos_erf(x):
    a1 =  0.254829592
    a2 = -0.284496736
    a3 =  1.421413741
    a4 = -1.453152027
    a5 =  1.061405429
    p  =  0.3275911
    t = 1.0 / (1.0 + p * x)
    poly = a1*t + a2*(t**2) + a3*(t**3) + a4*(t**4) + a5*(t**5)
    y = 1.0 - poly * np.exp(-x*x)
    return y

def golden(a, b, tol=1e-8):
    r = (np.sqrt(5)-1)/2    # Golden ratio
    x1 = a+(1-r)*(b-a)      # Compute left interior point
    x2 = a+r*(b-a)          # Compute right interior point
    f1 = f(x1)              # Function at x1
    f2 = f(x2)              # Function at x2
    while abs(b-a) > tol:           # Loop until tolerance is reached
        if f1 < f2:                 # If f1 is less than f2, the minimum is in [a, x2]
            b = x2                  
            x2 = x1                 
            f2 = f1                 # Update f2 to f1, since x1 moved to x2
            x1 = a+(1-r)*(b-a)      # Compute a new x1 in new interval [a, b]
            f1 = f(x1)              # Function at new x1
        else:                       # Otherwise, the minimum is in [x1, b]
            a = x1                  
            x1 = x2                 
            f1 = f2                 # Update f1 to f2, since x2 moved to x1
            x2 = a+r*(b-a)          # Compute a new x2 in the updated interval [a, b]
            f2 = f(x2)              # Function at new x2
    return a

if __name__ == "__main__":
    print(f"Approximate function minimum: {golden(0,2)}")

# Initial plot
# x = [i*0.02 for i in range(100)]
# y = [f(i) for i in x]
# plt.scatter(x, y)
# plt.grid()
# plt.show()