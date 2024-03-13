# --------------------- Numerical Root Solver ------------------- #

'''
This program allows the user to solve quadratic equations by numerical methods.
The function is plotted to allow the user to make a better guess of what the root will be.
'''

import matplotlib.pyplot as plt
import scipy.optimize as sci
import sympy as smp
import numpy as np

def solver(xmin, xmax, ymin, ymax, func, guess):

    xlist = np.linspace(xmin, xmax, 1000)

    plt.plot(xlist, func(xlist), color= 'cyan')
    plt.plot(xlist, xlist*0, linestyle= '--', color= 'red')
    plt.ylim(ymin, ymax)
    plt.show()

    roots = sci.fsolve(func, [guess])
    
    return roots
