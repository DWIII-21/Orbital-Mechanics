#---------- Converter for Cartesian Element to Classical Elements ----------#

'''
INPUTS
-------
Gravitation Parameter (mu)
rcart = [x,y,z] 
vcart = [xdot,ydot,zdot]
-------
OUTPUTS
------
Classical elements:
Semi-Major Axis (a)
Eccentricity (e)
Inclination (i)
Right Angle of the Ascending Node (omega)
Argument of Periapsis (w)
Mean Anomaly(M)
'''

import math as m
import numpy as np

def cart2class(mu, x, y, z, xdot, ydot, zdot):
    rcart = np.array([x,y,z])
    vcart = np.array([xdot,ydot,zdot])

    r = np.linalg.norm(rcart)
    v = np.linalg.norm(vcart)
    
    hcart = np.array(np.cross(rcart,vcart,axis=0))
    h = np.linalg.norm(hcart)
    xhat = np.array([1,0,0])
    zhat = np.array([0,0,1])  #Unit vector z
    zhat = np.reshape(zhat,(3,1))
    ahat = np.cross(zhat,hcart,axis=0)/np.linalg.norm(np.cross(zhat,hcart,axis=0)) #Unit vector a

    a = 1/(2/r-v**2/mu)
    ecart = (np.cross(vcart,hcart,axis=0))/mu - rcart/r
    e = np.linalg.norm(ecart)
    i = np.arccos(np.vdot(hcart,zhat)/h)
    omega = np.arccos(np.vdot(xhat,ahat))
    if ahat[1] >= 0:
        omega = omega
    elif ahat[1] < 0:
        omega = 2*m.pi-omega
    w = np.arccos(np.vdot(ahat,ecart)/e)
    if ecart[2] >= 0:
        w = w
    elif ecart[2] < 0:
        w = 2*m.pi-w
    nu = np.arccos(np.vdot(ecart,rcart)/(e*r))
    if np.vdot(rcart,vcart) >= 0:
        nu = nu
    elif np.vdot(rcart,vcart) < 0:
        nu = 2*m.pi-nu
    
    #E = 2*m.atan(m.tan(nu/2)/m.sqrt((1+e)/(1-e)))
    E = m.acos((a-r)/(a*e))
    if np.vdot(rcart,vcart) >= 0:
        E = E

    else:
        E = 2*m.pi-E
    M = E-e*m.sin(E)
    return a, e, i, omega, w, M
