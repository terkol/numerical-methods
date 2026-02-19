import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import BSpline

if __name__ == "__main__":
    nodes = [0, 1, 2, 3]    # Set nodes

    # Create B-Spline for ranges of nodes
    b_spline1 = BSpline.basis_element(nodes[:2])    
    b_spline2 = BSpline.basis_element(nodes[:3])
    b_spline3 = BSpline.basis_element(nodes[:4])

    # Set x axis for plot
    x = np.linspace(0, 3, 100)

    #Plot
    plt.plot(x, b_spline1(x), label="k=1")
    plt.plot(x, b_spline2(x), label="k=2")
    plt.plot(x, b_spline3(x), label="k=3")
    plt.legend()
    plt.grid()
    plt.show()

