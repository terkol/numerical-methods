import numpy as np
import os

def design_matrix(t):
    A = []
    for i in t:
        A.append([1, i, np.cos(i)])     # Assemble matrix row by row
    return np.array(A)

def solve_for_x(A, y):
    return np.linalg.solve(A.T @ A, A.T @ y)    # Use numpy's linear algebra solve function to find x

def moore_penrose(A, y):
    return np.linalg.pinv(A) @ y    # Use numpy's pseudoinverse function to find the pseudoinverse of A

def deviation(t, y, x):
    errors = []
    for i in range(len(t)):                                     # For every line in linear_data.txt
        errors.append(y[i]-(x[0]+x[1]*t[i]+x[2]*np.cos(t[i])))  # Calculate the deviation from y
    errors = np.array(errors)   # Change into a numpy array
    return np.sqrt(sum(errors**2)/(len(errors)-len(x)))         # Return standard deviation

def x_errors(A, sd, x):
    ATA_inv = np.linalg.inv(A.T @ A) # Invert A.T @ A for error propagation calculations
    errors = []
    errors.append(sd*np.sqrt((ATA_inv)[0][0]))  # Calculate error for α by multiplying inverted matrixes' index by the variance and take the square root
    errors.append(sd*np.sqrt((ATA_inv)[1][1]))
    gamma_square_error = sd*np.sqrt((ATA_inv)[2][2]) # Since gamma is squared...
    errors.append(gamma_square_error/(2*np.sqrt(x[2]))) # ...The error needs to also be propagated through the squareroot
    return errors

if __name__ == "__main__":
    data = np.loadtxt(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))+"\\data\\linear_data.txt")
    t_list = data[:,0]
    y_list = data[:,1]

    A = design_matrix(t_list)
    param = solve_for_x(A, y_list)
    print(f"Least squares: α={param[0]}, β={param[1]}, γ={np.sqrt(param[2])}")

    param = moore_penrose(A, y_list)
    print(f"Moore-Penrose: α={param[0]}, β={param[1]}, γ={np.sqrt(param[2])}")

    sd = deviation(t_list, y_list, param)
    print(f"Standard deviation of noise: {sd}")
    param_errors = x_errors(A, sd, param)
    print(f"α≈{str(param[0])[:6]}±{str(param_errors[0])[:6]}, "
        f"β≈{str(param[1])[:6]}±{str(param_errors[1])[:6]}, "
        f"γ≈{str(param[2])[:6]}±{str(param_errors[2])[:6]}")


