# SQP test
# Author: Spencer Boone
import numpy as np
import copy
import matplotlib.pyplot as plt
from scipy.optimize import minimize

def obj_function(x):
    obj = x**4 - x**3 - 6*x**2 + 4*x + 12
    return obj

def obj_gradient(x):
    grad_obj = 4*x**3 - 3*x**2 - 12*x + 4
    return grad_obj

def obj_hessian(x):
    hess_obj = 12*x**2 - 6*x - 12
    return hess_obj


# Initialize algorithm
x_0 = -1
tol = 1e-6
obj_init = obj_function(x_0)
obj_diff = np.abs(obj_init)
x = copy.deepcopy(x_0)
i = 0

# Initial plot
plt.figure()
x_plot = np.linspace(-3,3,1001)
plt.plot(x_plot, obj_function(x_plot), label = '$f(x)$')
plt.scatter(x, obj_init, 50, marker = '^', color = 'crimson', zorder = 3, label = 'Current guess for min[$f(x)$]')
plt.legend()
plt.show()
plt.ylim([-10,40])


# Run SQP
while obj_diff > tol:
    
    # Compute f(x), f', f''
    obj_iter = obj_function(x)
    obj_grad = obj_gradient(x)
    obj_hess = obj_hessian(x)
    
    # Compute dx
    dx = -obj_grad/obj_hess
    
    # Compute next iteration and compare objective function values
    x_new = x + dx
    obj_iter_new = obj_function(x_new)    
    obj_diff = np.abs(obj_iter_new - obj_iter)
    
    
    # Plot iteration
    plt.figure()
    plt.plot(x_plot, obj_function(x_plot), label = '$f(x)$')
    plt.plot(x_plot, obj_iter + obj_grad*(x_plot-x) + 1/2*obj_hess*(x_plot-x)**2, label = '2nd-order expansion')
    plt.scatter(x, obj_iter, 50, marker = '^', color = 'crimson', zorder = 3, label = 'Current guess for min[$f(x)$]')
    plt.scatter(x_new, obj_iter_new, 50, marker = '^', color = 'tab:olive', zorder = 3, label = 'Next guess for min[$f(x)$]')
    plt.scatter(x_new, obj_iter + obj_grad*(x_new-x) + 1/2*obj_hess*(x_new-x)**2, 50, marker = '^',
                color = 'tab:pink', zorder = 3, label = 'Minimum of 2nd-order expansion')
    plt.plot([x_new, x_new], [obj_iter_new, obj_iter + obj_grad*(x_new-x) + 1/2*obj_hess*(x_new-x)**2], 'r--')

    plt.show()
    plt.ylim([-10,60])
    plt.legend()
    x = copy.deepcopy(x_new)
    i = i + 1

    
 
# Test out python's SQP solver
x_scipy_opt = minimize(obj_function, x_0, method='SLSQP', tol=tol,
                       options = {'maxiter': 100})

# Test out python's SQP solver with gradient
x_scipy_opt_gradient = minimize(obj_function, x_0, method='SLSQP', tol=tol,
                                jac = obj_gradient, options = {'maxiter': 100})
