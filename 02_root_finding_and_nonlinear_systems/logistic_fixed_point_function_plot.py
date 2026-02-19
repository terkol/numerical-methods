import numpy as np
import matplotlib.pyplot as plt

def f(x, mu):
    return x - mu * x * (1 - x)

if __name__ == "__main__":
    mu_values = [2.5, 3.0, 3.2, 3.5]
    x_values = np.linspace(0.3, 1, 400) 


    for mu in mu_values:
        y_values = f(x_values, mu)
        plt.plot(x_values, y_values, label=f'μ = {mu}')

    plt.xlabel("x")
    plt.ylabel("f(x)")
    plt.legend()
    plt.grid(True)
    plt.show()
