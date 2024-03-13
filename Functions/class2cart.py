 #---------- Converter for Classical Elements to Cartesian Elements ----------#

'''
INPUTS
------
Gravitation Parameter (mu)
Classical elements:
Semi-Major Axis (a)
Eccentricity (e)
Inclination (i)
Right Angle of the Ascending Node (omega)
Argument of Periapsis (w)
Mean Anomaly (M)
-------
OUTPUTS
-------
rcart = [x,y,z] 
vcart = [xdot,ydot,zdot]
'''
import math as m
import numpy as np
import sympy as smp
import NewtonSolver as ns

def class2cart(mu,a,e,i,omega,w,M):
    p = a*(1-e**2)
    h = m.sqrt(mu*p)

    x = smp.symbols('x')
    y = x-e*smp.sin(x)-M
    E = ns.newton(y,smp.pi)

    nu = (m.atan(m.sqrt((1+e)/(1-e))*m.tan(E/2)))*2+2*m.pi
    nu_deg = nu*180/m.pi
    r = p/(1+e*m.cos(nu))
    rpol = np.array([r,0,0])
    rpol = rpol.reshape(3,1)
    vr = mu*e*m.sin(nu)/h 
    vo = mu/h*(1+e*m.cos(nu))
    vpol = np.array([vr,vo,0])
    vpol = vpol.reshape(3,1)

    theta = nu+w
    DCM = [m.cos(omega)*m.cos(theta)-m.sin(omega)*m.sin(theta)*m.cos(i), -m.cos(omega)*m.sin(theta)-m.sin(omega)*m.cos(theta)*m.cos(i), m.sin(omega)*m.sin(i)  ,
    m.sin(omega)*m.cos(theta)+m.cos(omega)*m.sin(theta)*m.cos(i), -m.sin(theta)*m.sin(omega)+m.cos(omega)*m.cos(theta)*m.cos(i), -m.cos(omega)*m.sin(i)  ,
    m.sin(theta)*m.sin(i), m.cos(theta)*m.sin(i), m.cos(i)]
    DCM = np.array(DCM)
    DCM = DCM.reshape((3,3))

    def cartConv(dcm,polarvector):
        cartesianvector = dcm @ polarvector  #In python '@' is the operator for matrix multiplication
        return cartesianvector
    
    rcart = cartConv(DCM,rpol)
    vcart = cartConv(DCM,vpol)

    return rcart,vcart
