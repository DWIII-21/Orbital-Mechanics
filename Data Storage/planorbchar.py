#---------------- Celestial Body Characteristic Calculator ---------------#

#This script performs two functions
#The first function (planchar) returns characteristics of the chosen celestial body
#The second function (orbchar) returns characteristics of the chosen celestial body's orbit around the Sun

import math as m
import numpy as np

###### SETUP ######


rad = np.pi/180  #radian conversion from degrees

### The characteristics of the planets themselves ###

#Mercury
mu_mc = 2.203*10**4
eq_mc = 2440  #Equitorial radius in km
M_mc = 0.330*10**24 #Mass of Mercury in kg
clr_mc = 'bisque'

#Venus
mu_vs = 3.2486*10**5
eq_vs = 6051.9
M_vs = 4.87*10**24
clr_vs = 'goldenrod'

#Earth
mu_eth = 3.986*10**5
eq_eth = 6378.14
M_eth = 5.9722*10**24
clr_eth = 'blue'

#Moon
mu_mn = 4.903*10**3
eq_mn = 1737.5
M_mn = 0.073*10**24
clr_mn = 'grey'

#Mars
mu_mr = 4.2828*10**4
eq_mr = 3390
M_mr = 0.642*10**24
clr_mr = 'indianred'

#Jupiter
mu_jp = 1.2669*10**8
eq_jp = 69911  #Jupiter has an equitorial radius of 69911 km
M_jp = 1898*10**24
clr_jp = 'darksalmon'

#Saturn
mu_sn = 3.7931*10**7
eq_sn = 58232  #Saturn has an equitorial radius of 58232 km
M_sn = 568*10**24
clr_sn = 'wheat'

#Titan
mu_tn = 8978.14
eq_tn = 2574.7
M_tn = 1.3455*10**23
clr_tn = 'mediumaquamarine'

#Uranus
mu_ur = 5.7940*10**6
eq_ur = 25362  #Uranus has an equitorial radius of 25362 km
M_ur = 86.8*10**24
clr_ur = 'paleturquoise'

#Neptune
mu_np = 6.8351*10**6
eq_np = 24624  #Neptune has an equitorial radius of 24624 km
M_np = 102*10**24
clr_np = 'royalblue'

#Pluto
mu_pt = 8.724*10**2
eq_pt = 1195  #Pluto has an equitorial radius of 1195 km
M_pt = 0.0130*10**24
clr_pt = 'rosybrown'

#Sun
muS = 1.327*10**11
eq_S = 696300
M_S = 1.9891*10**30
clr_S = 'orange'


### Determining characteristics for each planet's orbit about the Sun ###


#Mercury
a_mc = 0.387098*1.496*10**8   #A.U. conversion to km
e_mc = 0.205638
rp_mc = a_mc*(1-e_mc)
ra_mc = a_mc*(1+e_mc)
b_mc = a_mc*m.sqrt(1-e_mc**2)
p_mc = a_mc*(1-e_mc**2)
h_mc = m.sqrt(p_mc*muS)
r_mc = p_mc/(1+e_mc*m.cos(275*rad))

#Venus
a_vs = 0.723328*1.496*10**8
e_vs = 0.006753
rp_vs = a_vs*(1-e_vs)
ra_vs = a_vs*(1+e_vs)
b_vs = a_vs*m.sqrt(1-e_vs**2)
p_vs = a_vs*(1-e_vs**2)
h_vs = m.sqrt(p_vs*muS)
r_vs = p_vs/(1+e_vs*m.cos(120*rad))

#Earth
a_eth = 1.496*10**8
e_eth = 0.016882
rp_eth = a_eth*(1-e_eth)
ra_eth = a_eth*(1+e_eth)
b_eth = a_eth*m.sqrt(1-e_eth**2)
p_eth = a_eth*(1-e_eth**2)
h_eth = m.sqrt(p_eth*muS)
r_eth = p_eth/(1+e_eth*m.cos(70*rad))

#Earth's Moon
a_mn = 3.844*10**5
e_mn = 0.05490
rp_mn = a_mn*(1-e_mn)
ra_mn = a_mn*(1+e_mn)
b_mn = a_mn*m.sqrt(1-e_mn**2)
p_mn = a_mn*(1-e_mn**2)
h_mn = m.sqrt(p_mn*mu_eth)
r_mn = p_mn/(1+e_mn*m.cos(32*rad))

#Mars
a_mr = 1.523706*1.496*10**8
e_mr = 0.09352
rp_mr = a_mr*(1-e_mr)
ra_mr = a_mr*(1+e_mr)
b_mr = a_mr*m.sqrt(1-e_mr**2)
p_mr = a_mr*(1-e_mr**2)
h_mr = m.sqrt(p_mr*muS)
r_mr = p_mr/(1+e_mr*m.cos(250*rad))

#Jupiter
a_jp = 5.20149*1.496*10**8
e_jp = 0.048984
rp_jp = a_jp*(1-e_jp)
ra_jp = a_jp*(1+e_jp)
b_jp = a_jp*m.sqrt(1-e_jp**2)
p_jp = a_jp*(1-e_jp**2)
h_jp = m.sqrt(p_jp*muS)
r_jp = p_jp/(1+e_jp*m.cos(60*rad))

#Saturn
a_sn = 9.54327*1.496*10**8
e_sn = 0.005456
rp_sn = a_sn*(1-e_sn)
ra_sn = a_sn*(1+e_sn)
b_sn = a_sn*m.sqrt(1-e_sn**2)
p_sn = a_sn*(1-e_sn**2)
h_sn = m.sqrt(p_sn*muS)
r_sn = p_sn/(1+e_sn*m.cos(350*rad))

#Titan
a_tn = 1.221865*10**6
e_tn = 0.0288
rp_tn = a_tn*(1-e_tn)
ra_tn = a_tn*(1+e_tn)
b_tn = a_tn*m.sqrt(1-e_tn**2)
p_tn = a_tn*(1-e_tn**2)
h_tn = m.sqrt(p_tn*mu_tn)
r_tn = p_tn/(1+e_tn*m.cos(0*rad))

#Uranus
a_ur = 19.17113*1.496*10**8
e_ur = 0.048636
rp_ur = a_ur*(1-e_ur)
ra_ur = a_ur*(1+e_ur)
b_ur = a_ur*m.sqrt(1-e_ur**2)
p_ur = a_ur*(1-e_ur**2)
h_ur = m.sqrt(p_ur*muS)
r_ur = p_ur/(1+e_ur*m.cos(65*rad))

#Neptune
a_np = 29.99375*1.496*10**8
e_np = 0.008892
rp_np = a_np*(1-e_np)
ra_np = a_np*(1+e_np)
b_np = a_np*m.sqrt(1-e_np**2)
p_np = a_np*(1-e_np**2)
h_np = m.sqrt(p_np*muS)
r_np = p_np/(1+e_np*m.cos(0))

#Pluto
a_pt = 39.2305*1.496*10**8
e_pt = 0.247975
rp_pt = a_pt*(1-e_pt)
ra_pt = a_pt*(1+e_pt)
b_pt = a_pt*m.sqrt(1-e_pt**2)
p_pt = a_pt*(1-e_pt**2)
h_pt = m.sqrt(p_pt*muS)
r_pt = p_pt/(1+e_pt*m.cos(305*rad))


###### CALCULATIONS ######


# planchar Function: Returns characteristics of each planet/body
def planchar(planet):

    if planet == 'Mercury' or planet == 'mercury':
        mu = mu_mc
        eq = eq_mc
        M = M_mc
        clr = clr_mc

    elif planet == 'Venus' or planet == 'venus':
        mu = mu_vs
        eq = eq_vs
        M = M_vs
        clr = clr_vs

    elif planet == 'Earth' or planet ==  'earth':
        mu = mu_eth
        eq = eq_eth
        M = M_eth
        clr = clr_eth

    elif planet == 'Moon' or planet == 'moon':
        mu = mu_mn
        eq = eq_mn
        M = M_mn
        clr = clr_mn

    elif planet == 'Mars' or planet == 'mars':
        mu = mu_mr
        eq = eq_mr
        M = M_mr
        clr = clr_mr

    elif planet == 'Jupiter':
        mu = mu_jp
        eq = eq_jp
        M = M_jp
        clr = clr_jp

    elif planet == 'Saturn':
        mu = mu_sn
        eq = eq_sn
        M = M_sn
        clr = clr_sn

    elif planet == 'Titan':

        mu = mu_tn
        eq = eq_tn
        M = M_tn
        clr = clr_tn

    elif planet == 'Uranus':
        mu = mu_ur
        eq = eq_ur
        M = M_ur
        clr = clr_ur

    elif planet == 'Neptune':
        mu = mu_np
        eq = eq_np
        M = M_np
        clr = clr_np

    elif planet == 'Pluto':
        mu = mu_pt
        eq = eq_pt
        M = M_pt
        clr = clr_pt

    elif planet == 'Sun' or planet == 'sun':
        mu = muS
        eq = eq_S
        M = M_S
        clr = clr_S

    else:
        print("Please restart and choose from available directory")
        exit(1)

    return mu,eq,M,clr

# orbchar Function: Returns characteristics of each planet's orbit about the Sun
def planmisc(planet):
    if planet == 'Mercury' or planet == 'mercury':
        temp = 333 # Fahrenheit
        moons = 0
        gravity = 3.7 # in m/s^2
        inclination = 7.0 # in degrees
        day_length = 4222.6 # in Earth hours

    elif planet == 'Venus' or planet == 'venus':
        temp = 867
        moons = 0
        gravity = 8.9
        inclination = 3.4
        day_length = 2802.0

    elif planet == 'Earth' or planet ==  'earth':
        temp = 59
        moons = 1
        gravity = 9.8
        inclination = 0.0
        day_length = 24.0

    elif planet == 'Moon' or planet == 'moon':
        moons = 0
        gravity = 1.6
        inclination = 5.1
        day_length = 708.7

    elif planet == 'Mars' or planet == 'mars':
        temp = -85
        moons = 2
        gravity = 3.7
        inclination = 1.8
        day_length = 24.7

    elif planet == 'Jupiter':
        temp = -166
        moons = 79
        gravity = 23.1
        inclination = 1.3
        day_length = 9.9

    elif planet == 'Saturn':
        temp = -220
        moons = 82
        gravity = 9.0
        inclination = 2.5
        day_length = 10.7

    elif planet == 'Titan':
        pass

    elif planet == 'Uranus':
        temp = -320
        moons = 27
        gravity = 8.7
        inclination = 0.8
        day_length = 17.2

    elif planet == 'Neptune':
        temp = -330
        moons = 14
        gravity = 11.0
        inclination = 1.8
        day_length = 16.1

    elif planet == 'Pluto':
        temp = -375
        moons = 5
        gravity = 1.3
        inclination = 17.2
        day_length = 153.3

    elif planet == 'Sun' or planet == 'sun':
        temp = 10000
        moons = 0
        gravity = 'n/a'
        inclination = 'n/a'
        day_length = 'n/a'

    else:
        print("Please restart and choose from available directory")
        exit(1)

    return temp,moons,gravity,inclination,day_length

# orbchar Function: Returns characteristics of each planet's orbit about the Sun
def orbchar(planet):

    if planet == 'Mercury' or planet == 'mercury':
        a = a_mc
        e = e_mc
        rp = rp_mc
        ra = ra_mc
        b = b_mc
        p = p_mc
        h = h_mc

    elif planet == 'Venus' or planet == 'venus':
        a = a_vs
        e = e_vs
        rp = rp_vs
        ra = ra_vs
        b = b_vs
        p = p_vs
        h = h_vs

    elif planet == 'Earth' or planet == 'earth':
        a = a_eth
        e = e_eth
        rp = rp_eth
        ra = ra_eth
        b = b_eth
        p = p_eth
        h = h_eth

    elif planet == 'Moon' or planet == 'moon':
        a = a_mn
        e = e_mn
        rp = rp_mn
        ra = ra_mn
        b = b_mn
        p = p_mn
        h = h_mn

    elif planet == 'Mars' or planet == 'mars':
        a = a_mr
        e = e_mr
        rp = rp_mr
        ra = ra_mr
        b = b_mr
        p = p_mr
        h = h_mr

    elif planet == 'Jupiter':
        a = a_jp
        e = e_jp
        rp = rp_jp
        ra = ra_jp
        b = b_jp
        p = p_jp
        h = h_jp

    elif planet == 'Saturn':
        a = a_sn
        e = e_sn
        rp = rp_sn
        ra = ra_sn
        b = b_sn
        p = p_sn
        h = h_sn

    elif planet == 'Titan':
        a = a_tn
        e = e_tn
        rp = rp_tn
        ra = ra_tn
        b = b_tn
        p = p_tn
        h = h_tn

    elif planet == 'Uranus':
        a = a_ur
        e = e_ur
        rp = rp_ur
        ra = ra_ur
        b = b_ur
        p = p_ur
        h = h_ur

    elif planet == 'Neptune':
        a = a_np
        e = e_np
        rp = rp_np
        ra = ra_np
        b = b_np
        p = p_np
        h = h_np

    elif planet == 'Pluto':
        a = a_pt
        e = e_pt
        rp = rp_pt
        ra = ra_pt
        b = b_pt
        p = p_pt
        h = h_pt

    else:
        print("Please restart and choose from available directory")
        exit(1)

    return a,e,rp,ra,b,p,h
