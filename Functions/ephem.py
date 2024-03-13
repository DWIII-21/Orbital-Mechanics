#---------- JPL Ephemeris Data Function ----------#

'''
INPUTS
-------
JDi
Body
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
import sympy as smp
import jplephem as jpl
from jplephem.spk import SPK
import NewtonSolver as ns
import cart2class as ca


def ephem(JDi,BODY):
    kernel = SPK.open('de421.bsp')
    np.set_printoptions(precision=3)

    AUtoKM = 1.49597870691*10**8
    #JDS = 2456594
    JDS = 2414864.50
    #JDM = 2463244
    JDM = 2471184.50
    body = 0
    cont = 0
    tol = 1*10**(-15)
    mu = 1.32712440018*10**11  #Gravitational Parameter of the Sun

    if (JDS < JDi) and (JDi <=JDM):
        cont = 1
        JDe = m.floor(JDi - JDS); #Whole number component
        JDd = JDi - m.floor(JDi); #Decimal component
    else:
        print('Julian Date not in accepted range')
        return 
    
    dT = JDd*86400

    if BODY == 'Mercury':
        bary = 1
        body = 199
    elif BODY == 'Venus':
        bary = 2
        body = 299
    elif BODY == 'Earth':
        bary = 3
        body = 399
    elif BODY == 'Moon':
        bary = 3
        body = 301
    elif BODY == 'Mars':
        bary = 4
        body = 499

    position1, velocity1 = kernel[0,bary].compute_and_differentiate(JDi)  #Finds position and velo of BODY's barycenter around the Sun
    position2, velocity2 = kernel[bary,body].compute_and_differentiate(JDi) #Finds poisiton and velo of BODY w/rt it's barycenter

    position = position1 - position2  #Calculates velo of BODY w/rt the Sun in KM
    velocity = velocity1 - velocity2  #Calculates velo of BODY w/rt the Sun in KM/DAY
    velocity_per_second = velocity/86400  #Translates velo into KM/SECOND
    rcart = position
    vcart = velocity_per_second

    r = np.linalg.norm(rcart)
    v = np.linalg.norm(vcart)
    
    hcart = np.array(np.cross(rcart,vcart,axis=0))
    h = np.linalg.norm(hcart)
    xhat = np.array([1,0,0])
    zhat = np.array([0,0,1])  #Unit vector z
    zhat = np.reshape(zhat,(3,1))
    ahat = np.cross(zhat,hcart,axis=0)/np.linalg.norm(np.cross(zhat,hcart,axis=0)) #Unit vector a

    a = 1/(2/r-v**2/mu)
    n = m.sqrt(mu/a**3)
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
    
    E = m.acos((a-r)/(a*e))
    if np.vdot(rcart,vcart) >= 0:
        E = E
    else:
        E = 2*m.pi-E
    M = E-e*m.sin(E)

    Mf = M + dT*n
    Eg = Mf + e*m.sin(Mf) + (e**2/2)*m.sin(2*Mf)

    if 0 < JDd:
        dif = Eg - Mf
        k = 1
        kmax = 15
        x = smp.symbols('x')
        fx = x-e*smp.sin(x)
        kn = ns.newton(fx,k)
    
        Ef = Eg
        dE = Ef-E
        rf = a*(1-e*m.cos(Ef))

        f = 1 - (a/r)*(1-m.cos(dE))
        g = dT - m.sqrt(a**3/mu)*(dE-m.sin(dE))
        fdot = (-m.sqrt(mu*a)*m.sin(dE))/(rf*r)
        gdot = 1 - (a/rf)*(1-m.cos(dE))
    
        rfinal = f*rcart + g*vcart
        vfinal = fdot*rcart + gdot*vcart
    else:
        rfinal = rcart
        vfinal = vcart
    
    r = np.linalg.norm(rfinal)
    v = np.linalg.norm(vfinal)

    if np.vdot(rfinal,vfinal) >= 0:
        E = m.acos((a-r)/(a*e))
    else:
        E = 2*m.pi - m.acos((a-r)/(a*e))  

    M = E - e*m.sin(E)

    kernel.close()
    return a,e,i,omega,w,M
