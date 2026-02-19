import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline

if __name__ == "__main__":
    # Generate random data with constant seeds and a set of nodes for each
    np.random.seed(1)
    y5 = np.random.rand(5)
    x5 = np.array(range(5))*4
    np.random.seed(2)
    y10 = np.random.rand(10)
    x10 = np.array(range(10))*2
    np.random.seed(3)
    y20 = np.random.rand(20)
    x20 = list(range(20))

    # Create splining instances with different boundary conditions
    spline_bc0 = CubicSpline(x5, y5, bc_type=((2, 0), (2, 0)))
    spline_bc1 = CubicSpline(x10, y10, bc_type=((2, 1.0), (2, 1.0)))
    spline_bc10 = CubicSpline(x20, y20, bc_type=((2, 10.0), (2, 10.0)))

    # Spline random data
    x5_spline = np.linspace(x5[0], x5[-1], 200)
    y5_spline = spline_bc0(x5_spline)
    x10_spline = np.linspace(x10[0], x10[-1], 200)
    y10_spline = spline_bc1(x10_spline)
    x20_spline = np.linspace(x20[0], x20[-1], 200)
    y20_spline = spline_bc10(x20_spline)

    # Plot splined data
    plt.scatter(x5, y5,label=5,color='g')
    plt.plot(x5_spline, y5_spline,color='g')
    plt.scatter(x10, y10,label=10,color='r')
    plt.plot(x10_spline, y10_spline,color='r')
    plt.scatter(x20, y20,label=20,color='b')
    plt.plot(x20_spline, y20_spline,color='b')
    plt.legend()
    plt.show()