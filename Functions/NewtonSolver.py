 #---------- Newton's Method Solver ----------#

import math as m
import numpy as np
import sympy as smp
from sympy import *

def newton(f,x0):
    '''Approximate solution of f(x)=0 by Newton's method.

    Parameters
    ----------
    f : function
        Function for which we are searching for a solution f(x)=0.
    df : function
        Derivative of f(x).
    x0 : number
        Initial guess for a solution f(x)=0.
    tolerance : number
        Stopping criteria is abs(f(x)) < tolerance.
        Currently set to 10^-6
    max_iter : integer
        Maximum number of iterations of Newton's method.
        Currently set to 100

    Returns
    -------
    xn : number
        Implement Newton's method: compute the linear approximation
        of f(x) at xn and find x intercept by the formula
            x = xn - f(xn)/df(xn)
        Continue until abs(f(xn)) < tolerance and return xn.
        If df(xn) == 0, return None. If the number of iterations
        exceeds max_iter, then return None.

    Example (SOLVING KEPLER'S EQUATION)
    --------
    >>> x = smp.symbols('x')  (IMPORTANT)
    >>> y = x - e*smp.sin(x) - M  (General)
    >>> y = x - 0.23*smp.sin(x) - 250*smp.pi/180 (Example)
    >>> E = newton(y,smp.pi)
    Found solution after 3 iterations.
    4.166723272394594  (measured in Radians)
    '''

    df = smp.diff(f)
    tolerance = 10**(-6)
    max_iter = 100

    x = smp.symbols('x')
    xn = x0
    for n in range(0,max_iter):
        fxn = f.subs(x,xn)
        if abs(fxn) < tolerance:
            # print('Found solution after',n,'iterations.')
            return float(xn)
        dfxn = df.subs(x,xn)
        if dfxn == 0:
            # print('Zero derivative. No solution found.')
            return None
        xn = xn - fxn/dfxn
    # print('Exceeded maximum iterations. No solution found.')
    return None
