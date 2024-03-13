# ------------------------ GEOCENTRIC ORBIT PROP TESTER --------------------------- #

import numpy as np
import planetary_data as pd

import geocentric_orbit_propagator as geop


## TESTING GEOCENTRIC ORBIT PROP ##

def main():

    geop.plot_orbits([[pd.earth['eq'] + 417, 0.0001788, np.radians(51.6398), np.radians(207.9798), np.radians(250.1293), np.radians(148.6385), 'classical']],
                dt= 5, tspan= 1e4, labels= ['ISS Orbit'], perifocal= True, trackm= 120) 

    geop.plot_orbits([[26528.14, 0.7369532880933227, 1.106538745764405, 0.5235987755982988, 4.71238898038469, 3.141592653589793, 'classical']],
                dt= 5, tspan= 7.5e4, labels= ['Molniya Orbit'], perifocal= True, trackh= 12)

    geop.plot_orbits([[26560, 0.01, np.radians(55), 0, np.radians(30), 0, 'classical']],
                dt= 5, tspan= 1e5, labels= ['GPS Satellite Orbit'], perifocal= True, trackh= 4.5)
    
if __name__ == '__main__':
    main()
