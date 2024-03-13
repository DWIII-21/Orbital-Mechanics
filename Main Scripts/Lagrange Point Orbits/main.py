# --------------------- Lagrange Orbits w/ CR3BP ----------------------- #

'''
Plots orbits of Earth-Moon Lagrange points to determine their stability.
'''

import numpy as np
import Orbital_Propagator as prop
import planetary_data as pd
import matplotlib.pyplot as plt
import scipy.optimize as sci
import sympy as smp
import numerical_root_solver as rs

def main():

    ####################
    # Earth-Moon System #
    ####################

    dt = 0.005
    tspan = 30
    dx = -0.01  # x-position pertabation to determine stability of each Lagrange Point (only L4 and L5 are stable)

    names = ['L1', 'L2', 'L3', 'L4', 'L5']
    colors = ['magenta', 'steelblue', 'lime', 'cyan', 'crimson']

    EM = prop.CR3BP_system('earth-moon')
    mu = EM['mu']

    lagrange_orbit1 = prop.CR3BP(system= 'em', dt= dt, tspan= tspan, y0= [EM['L1'], 0, 0, 0, 0, 0], dx= dx, objlabel= names[0], maxstep= 0.05)
    lagrange_orbit2 = prop.CR3BP(system= 'em', dt= dt, tspan= tspan, y0= [EM['L2'], 0, 0, 0, 0, 0], dx= dx, objlabel= names[1], maxstep= 0.05)
    lagrange_orbit3 = prop.CR3BP(system= 'em', dt= dt, tspan= 45, y0= [EM['L3'], 0, 0, 0, 0, 0], dx= dx, objlabel= names[2], maxstep= 0.02)
    lagrange_orbit4 = prop.CR3BP(system= 'em', dt= dt, tspan= tspan, y0 = [1/2 - mu, np.sqrt(3)/2, 0, 0, 0, 0], dx= dx, objlabel= names[3])
    lagrange_orbit5 = prop.CR3BP(system= 'em', dt= dt, tspan= tspan, y0 = [1/2 - mu, -np.sqrt(3)/2, 0, 0, 0, 0], dx= dx, objlabel= names[4])

    tlist = []
    ylist = []
    orbits = [lagrange_orbit1, lagrange_orbit2, lagrange_orbit3, lagrange_orbit4, lagrange_orbit5]

    # Plotting results for CR3BP
    plt.style.use('dark_background')
    plt.figure()
    plt.plot(-mu, 0, marker= '.', markersize= 18, color= lagrange_orbit1.b1['clr'], label= lagrange_orbit1.b1['name'])     # Plotting Large Body (Body 1)
    plt.plot(1 - mu, 0, marker= '.', markersize= 15, color= lagrange_orbit1.b2['clr'], label= lagrange_orbit1.b2['name'])   # Plotting Smaller Body (Body 2)

    n = 0
    for orbit in orbits:

        t, y = orbit.propagate_orbit()
        tlist.append(t)
        ylist.append(y)

        plt.plot( y[0, 0], y[0, 1], marker= '.', markersize= 10, color= colors[n], label= f'{names[n]}')
        plt.plot(y[:, 0], y[:, 1], color= colors[n], label= f'Trajectory of {names[n]}')

        n += 1

    plt.grid(linestyle= '--', linewidth= 0.2)
    plt.tight_layout()
    plt.title('Lagrange Point Stability for Earth-Moon System')
    plt.legend(bbox_to_anchor=(0.95, 0.7), loc='lower center')
    plt.show()


if __name__ == '__main__':
    main()
