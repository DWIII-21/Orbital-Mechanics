#----------------- C3 Calculator ----------------#

import numpy as np
import ephem as eph
import lambert as lmb
import class2cart as cl
import planorbchar as po
import planetary_data as pd

'''
INPUTS
-------
JDdp
JDar
Gravitational Parameter (mu)
Body 1 (as a string)
Body 2 (as a string)
Altitude above Body 2 at which incoming velocity is measured
-------
OUTPUTS
------
Characteristic Energy (C3) required to leave Body 1
Arrival velocity into Body 2's SOI at a desired altitude "h" (vinfm)
'''

def C3calc(JDdp,JDar,mu,body1,body2,h):

    [mu2,req2,mass2,color2] = po.planchar(body2)  # Finding gravitational parameter and equitorial radius of Body 2

    [a1,e1,i1,omega1,w1,M1] = eph.ephem(JDdp,body1)
    r1,v1 = cl.class2cart(mu,a1,e1,i1,omega1,w1,M1)
    [a2,e2,i2,omega2,w2,M2] = eph.ephem(JDar,body2)
    r2,v2 = cl.class2cart(mu,a2,e2,i2,omega2,w2,M2)

    # Rearranging v1 and v2 to allow proper vector subtraction
    v1 = v1.reshape((1,3))
    v2 = v2.reshape((1,3))

    [a,p,vDP,vAR,conv] = lmb.lambert(r1, r2, (JDar-JDdp)*24*60*60, mu, 1, 100, 10**(-14), 100)
    vDP = vDP.reshape((1,3))
    vAR = vAR.reshape((1,3))

    # Computing incoming and outgoing V-infinity values
    vinfp = np.linalg.norm(vDP - v1)
    vinfm = np.linalg.norm(vAR - v2)
    v_h = np.sqrt(2*(vinfm**2/2 + mu2/(req2+h))) # Incoming velocity at the desired altitude

    # Computing C3
    C3 = vinfp**2

    return C3,v_h
