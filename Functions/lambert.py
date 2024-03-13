#-------------- Lambert's Blackbox ---------------#
#
#   Programmer: Grant Hecht
#   Date:       3/18/2019
#   File:       lambert.py
#   Purpose:    This file contains the function lambert() which
#               solves Lambert's Problem using Battin's Method.
#

import numpy as np
import math

#
#   INPUTS
#   ------
#   R1:     Numpy Array of position vector of Pt#1 (in km)
#   R2:     Numpy Array of position vector of Pt#2 (in km)
#   TOF:    Time of Flight [sec]
#   mu:     Gravitational parameter for central body
#   JJ:     Integer that determines initial guess for x
#               (set JJ = 1 for an ellipse)
#               (set JJ = 0 for a parabola or hyperbola)
#   n:      Number of continued fraction levels
#               (set n = 0 for default of 10)
#               (recomend setting n >= 100 for most accurate results)
#   tol:    Tolerance to exit iterations
#                (set tol = 10^-14)
#   kmax:   Maximum number of iterations
#               (set kmax = 1000)
#   -------
#   OUTPUTS
#   -------
#   A:      Semi-Major Axis
#   P:      Semi-Latus Rectum
#   V1:     Numpy Array of Velocity Vector at Pt#1
#   V2:     Numpy Array of Velocity Vector at Pt#2
#   conv:   Boolean to indicate convergence
#
def lambert(R1,R2,TOF,mu,JJ,n,tol,kmax):

    # If does not converge,set as false
    conv = True

    # Sets n to default value if n = 0 is passed
    if n == 0:
        n = 10

    # Calculates Magnitude of R1 and R2
    r1 = np.linalg.norm(R1)
    r2 = np.linalg.norm(R2)

    # Find Transfer Angle
    ta = math.acos(np.vdot(R1,R2)/(r1*r2))
    if R1[0]*R2[1]-R1[1]*R2[0] < 0.:
        ta = 2*np.pi - ta

    # Find Chord
    c = math.sqrt(r1**2+r2**2-2.*r1*r2*math.cos(ta))

    # Find Semi-Perimeter
    s = (r1+r2+c)/2.

    # Find Lambda
    lamb = math.sqrt(r1*r2)*math.cos(ta/2.)/s

    # Find w
    w = math.atan(pow(r2/r1,0.25))-(np.pi/4.)

    # Find l
    if (0 < ta) and (ta < np.pi):
        l = (math.sin(ta/4.)**2+math.tan(2.*w)**2)/(math.sin(ta/4.)**2\
            +math.tan(2.*w)**2+math.cos(ta/2.))
    elif (np.pi <= ta) and (ta < 2.*np.pi):
        l = (math.cos(ta/4.)**2+math.tan(2.*w)**2-math.cos(ta/2.))/(\
            math.cos(ta/4.)**2+math.tan(2.*w)**2)
    else:
        print("Cannot Compute for Transfer Angle of 0 or 360 degrees")

    # Find m
    m = (8.*mu*TOF**2)/((s**3)*((1.+lamb)**6))

    # For Eliptical Transfer Orbit Use x = l for Initial Guess
    # For Hyperbolic or Parabolic Transfer Orbit Use x = 0
    if JJ == 0:
        x = 0.0
    else:
        x = l

    # Define Velocity Vectors
    V1 = np.array([0.0,0.0,0.0])
    V2 = np.array([0.0,0.0,0.0])
    y = 1

    # Define Delta x and Counter
    DX = 100.0
    k  = 0

    while abs(DX) > tol:
        # Breaks Itteration and sets conv = false if k = kmax
        if k >= kmax:
            conv = False
            break

        # Calculates Continued Fraction PHI for 'n' Levels
        eta = x/((math.sqrt(1.+x)+1.)**2)
        f = 1.0
        # Itterates to Calculate Levels 4 -> n
        for i in range(n,3,-1):
            ceta = (i**2)/((2.*i)**2-1.)
            f = 1.0 + ceta*eta/f
        # Finishes Calculations of PHI with levels 1 -> 3
        PHI = 8.*(math.sqrt(1.+x)+1.)/(3.+1./(5.+eta+(9.0/7.0)\
              *eta/f))

        h1 = ((l+x)**2)*(1.+3.*x+PHI)/((1.+2.*x+l)*(4.*x+PHI*(3.+x)))
        h2 = m*(x-l+PHI)/((1.+2.*x+l)*(4.*x+PHI*(3.+x)))
        B  = 27.*h2/(4.*(1.+h1)**3)
        u  = B/(2.*(math.sqrt(1.+B)+1.))

        # Calculates Continued Fraction K(u) for 'n' Levels
        f = 1.0
        r = int(n/2-1)
        # Itterate to Calculate Levels 3 -> n
        for j in range(r,0,-1):
            g2n  = 2.*(3.*j+1.)*(6.*j-1.)/(9.*(4.*j-1.)*(4.*j+1.))
            g2n1 = 2.*(3.*j+2.)*(6.*j+1.)/(9.*(4.*j+1.)*(4.*j+3.))
            f = 1. + g2n*u/(1. + g2n1*u/f)
        # Finishes Calculation of K(u) with levels 1 -> 2
        K = (1.0/3.0)/(1+(4.0/27.0)*u/f)

        # Calculates New Values for y and x
        yNew = ((1.+h1)/3.)*(2.+math.sqrt(1.+B)/(1.+2.*u*K**2))
        xNew = math.sqrt(((1.-l)/2)**2+m/(yNew**2))-(1.+l)/2.

        # Compared xNew with x
        DX = abs(xNew-x)

        # Sets x and y to xNew and yNew
        x = xNew
        y = yNew

        k = k+1

    # Computes Orbit Parameters and Initial and Final Velocity
    A = m*s*(1.+lamb)**2/(8.*x*y**2)
    P0 = (c**2)*(1.+x)**2/(16.*A*x)
    P = 4.*r1*r2*P0*math.sin(ta/2.)**2/c**2
    F = 1.-(r2/P)*(1.-math.cos(ta))
    G = r1*r2*math.sin(ta)/math.sqrt(mu*P)
    FDOT = math.sqrt(mu/P)*math.tan(ta/2.)*((1.-math.cos(ta))/P\
           -1./r2-1./r1)
    GDOT = 1.-(r1/P)*(1.-math.cos(ta))
    for i in range(0,3):
        V1[i]=(1/G)*(R2[i]-(F*R1[i]))
        V2[i]=FDOT*R1[i]+GDOT*V1[i]

    return (A,P,V1,V2,conv)    
