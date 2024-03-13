#---------------- Dictionary for Planetary Data ---------------#

'''
-------
VALUES LISTED
------
Name of Planet
Gravitational Parameter of Planet (km^3/s^2), Note: μ = M*G
Equitorial Radius of Planet (km)
Mass of Planet (kg)
Color of Planet
Colormap for Planet (for 3D plots)
Numeric ID - Body (used for gathering SPICE data of planet location; not applicable to all planets)
Numeric ID - Bary (used for gathering SPICE data of planet's barycenter location)
'''

#### To call values from this dictionary:
#       mu = planetary_data.sun['mu']

# Sun
sun = {
    'name': 'Sun',
    'mu': 1.327*10**11,
    'eq': 696300,
    'mass': 1.9891*10**30,
    'temp' : 10000,
    'moons' : 0,
    'grav' : 247,
    'inc' : 'n/a',
    'lod' : 'n/a',
    'clr': 'orange',
    'clrmap': 'hot',
    'id': 10
}

# Mercury
mercury = {
    'name': 'Mercury',
    'mu': 2.203*10**4,
    'eq': 2440,
    'mass': 0.330*10**24,
    'temp' : 333, # average temperature in Fahrenheit
    'moons' : 0,
    'grav' : 3.7, # gravity in m/s^2
    'inc' : 7.0, # inclination of orbit to ecliptic in degrees
    'lod' : 4222.6, # length of one day in Earth hours
    'a': 0.387098, # Semi-major axis in AU (1 AU = 1.496e8)
    'ecc': 0.205638,
    'period': 0.2408, # Orbital period of Mercury in years
    'clr': 'bisque',
    'clrmap': 'pink',
    'idbody': 199,
    'idbary': 1
}

# Venus
venus = {
    'name': 'Venus',
    'mu': 3.2486*10**5,
    'eq': 6051.9,
    'mass': 4.87*10**24,
    'temp' : 867,
    'moons' : 0,
    'grav' : 8.9,
    'inc' : 3.4, 
    'lod' : 2802.0,
    'a': 0.723328,
    'ecc': 0.006753,
    'period': 0.6152,
    'clr': 'goldenrod',
    'clrmap': 'autumn',
    'idbody': 299,
    'idbary': 2
}

# Earth
earth = {
    'name': 'Earth',
    'mu': 3.986*10**5,
    'eq': 6378.14,
    'mass': 5.9722*10**24,
    'temp' : 59,
    'moons' : 1,
    'grav' : 9.8,
    'inc' : 0.0,
    'lod' : 24.0,
    'a': 1,
    'ecc': 0.016882,
    'period': 1,
    'clr': 'blue',
    'clrmap': 'Blues',
    'idbody': 399,
    'idbary': 3
}

# Moon
moon = {
    'name': 'Moon',
    'mu': 4.903*10**3,
    'eq': 1737.5,
    'mass': 0.073*10**24,
    'moons' : 0,
    'grav' : 1.6,
    'inc' : 5.1,    # with respect to Earth's orbit
    'lod' : 708.7,
    'a': 3.844e5, # Semi-major axis around Earth IN KM!
    'ecc': 0.05490,
    'period': 27.322, # Orbital period in DAYS
    'clr': 'silver',
    'clrmap': 'Greys',
    'idbody': 301,
    'idbary': 3
}

# Mars
mars = {
    'name': 'Mars',
    'mu': 4.2828*10**4,
    'eq': 3390,
    'mass': 0.642*10**24,
    'temp' : -85,
    'moons' : 2,
    'grav' : 3.7,
    'inc' : 1.8,
    'lod' : 24.7,
    'a': 1.523706,
    'ecc': 0.09352,
    'period': 1.8808,
    'clr': 'crimson',
    'clrmap': 'Reds',
    'idbody': 499,
    'idbary': 4
}

# Jupiter
jupiter = {
    'name': 'Jupiter',
    'mu': 1.2669*10**8,
    'eq': 69911,
    'mass': 1898*10**24,
    'temp' : -166,
    'moons' : 79,
    'grav' : 23.1,
    'inc' : 1.3,
    'lod' : 9.9,
    'a': 5.20149,
    'ecc': 0.048984,
    'period': 11.86,
    'clr': 'darksalmon',
    'clrmap': 'copper'
}

# Europa
europa = {
    'name': 'Europa',
    #'mu': 1.2669*10**8,
    'eq': 1560.8,
    #'mass': 1898*10**24,
    #'temp' : -166,
    'moons' : 0,
    'grav' : 1.315,
    #'inc' : 1.3,
    #'lod' : 9.9,
    #'a': 5.20149,
    #'ecc': 0.048984,
    #'period': 11.86,
    'clr': 'papayawhip',
    #'clrmap': 'copper'
}

# Saturn
saturn = {
    'name': 'Saturn',
    'mu': 3.7931*10**7,
    'eq': 58232,
    'mass': 568*10**24,
    'temp' : -220,
    'moons' : 82,
    'grav' : 9.0,
    'inc' : 2.5,
    'lod' : 10.7,
    'a': 9.54327,
    'ecc': 0.005456,
    'period': 29.45,
    'clr': 'wheat',
    'clrmap': 'Wistia'
}

# Titan
titan = {
    'name': 'Titan',
    'mu': 8978.14,
    'eq': 2574.7,
    'mass': 1.3455*10**23,
    'clr': 'mediumaquamarine',
    'clrmap': 'summer'
}

# Uranus
uranus = {
    'name': 'Uranus',
    'mu': 5.7940*10**6,
    'eq': 25362,
    'mass': 86.8*10**24,
    'temp' : -320,
    'moons' : 27,
    'grav' : 8.7,
    'inc' : 0.8,
    'lod' : 17.2,
    'a': 19.17113,
    'ecc': 0.048636,
    'period': 84.07,
    'clr': 'paleturquoise',
    'clrmap': 'cool'
}

# Neptune
neptune = {
    'name': 'Neptune',
    'mu': 6.8351*10**6,
    'eq': 24624,
    'mass': 102*10**24,
    'temp' : -330,
    'moons' : 14,
    'grav' : 11.0,
    'inc' : 1.8,
    'lod' : 16.1,
    'a': 29.99375,
    'ecc': 0.008892,
    'period': 164.9,
    'clr': 'royalblue',
    'clrmap': 'winter'
}

# Pluto
pluto = {
    'name': 'Pluto',
    'mu': 8.724*10**2,
    'eq': 1195,
    'mass': 0.0130*10**24,
    'temp' : -375,
    'moons' : 5,
    'grav' : 1.3,
    'inc' : 17.2,
    'lod' : 153.3,
    'a': 39.2305,
    'ecc': 0.247975,
    'period': 249.6,
    'clr': 'rosybrown',
    'clrmap': 'pink'
}

# Spacecraft
spacecraft = {
    'name': 'Spacecraft',
    'clr': 'red',
    'clrmap': 'pink'
}
