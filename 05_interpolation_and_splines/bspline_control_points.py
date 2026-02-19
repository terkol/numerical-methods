import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import BSpline

if __name__ == "__main__":
    control_points = np.array([[0,0],[1,1],[2,0.5],[3,1],[4,0.5],[5,1.5],[6,1],[7,0.5],[8,0],[9,1.5],[10,0.5],[11,1],[12,1.5],[13,0]]) # Manually select 14 control points
    x = control_points[:, 0]    # Get x values of points
    y = control_points[:, 1]    # Get y values

    k = 3  # Choose degree
    knots = np.concatenate((np.zeros(k),np.linspace(0,1, len(control_points)-k+1), np.ones(k))) # Assemble knot vector

    spline_x = BSpline(knots, x, k) # Make spline for x
    spline_y = BSpline(knots, y, k) # Make spline for y

    x_plot = np.linspace(0, 1, 200) # Set frequency of x for forming plot
    t = spline_x(x_plot)    # Create spline for x axis
    C = spline_y(x_plot)    # Create spline for y axis

    # Plot
    plt.plot(t, C, label='Cubic B-spline')
    plt.plot(x, y, "o--", label='Control Points')
    plt.legend()
    plt.grid(True)
    plt.show()
