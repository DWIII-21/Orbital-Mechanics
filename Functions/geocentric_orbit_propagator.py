# ---------------------- Earth Orbit Visualizer --------------------- #

import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import *
from tkinter import ttk

import planetary_data as pd
import Orbital_Propagator as prop
import groundtracks as gt
import class2cart as cl
import cart2class as ct
import reference_frame_rotator as rfr



def plot_orbits(orbits, dt, tspan, labels= None, **kwargs):    # Inputs are the arrays containing orbit latlongs and the type of latlongs being input (classical or cartesian)

    disp_peri = False
    if 'perifocal' in kwargs:
        disp_peri = kwargs['perifocal']

    track = False
    if 'tracks' in kwargs:  # Looks for time in seconds
        track = True 
        ttt = kwargs['tracks']       # ttt = Time elapsed (given in seconds)
        time_conversion = 1
        time_units = 'Seconds'
    elif 'trackm' in kwargs:  # Looks for time in minutes
        track = True 
        ttt = kwargs['trackm'] * 60       # ttt = Time elapsed (given in minutes, converted on this line to seconds)
        time_conversion = 1 / 60
        time_units = 'Minutes'
    elif 'trackh' in kwargs:    # Looks for time in hours
        track = True 
        ttt = kwargs['trackh'] * 3600       # ttt = Time elapsed (given in hours, converted on this line to seconds)
        time_conversion = 1 / 3600
        time_units = 'Hours'
    elif 'trackd' in kwargs:    # Looks for time in days
        track = True 
        ttt = kwargs['trackd'] * 3600 * 24       # ttt = Time elapsed (given in days, converted on this line to seconds)
        time_conversion = 1 / 3600 / 24
        time_units = 'Days'

    save_fig = False
    if 'savefig' in kwargs:
        save_fig = True
        save_as = kwargs['savefig']

    steps = int(tspan/dt)
    n = 0
    for orbit in orbits:
        orbitcoords = orbit[:6]
        coord = orbit[6]
        if coord == 'cartesian':
            propped_orbit = prop.rk4_2Body(cb= 'earth', dt= dt, tspan = tspan, y0= orbitcoords, show= False)
            t, ys = propped_orbit.propagate_orbit()
            yin = ys[:, :3]

            a, e, inc, RAAN, w, M = ct.cart2class(pd.earth['mu'], orbitcoords[0], orbitcoords[1], orbitcoords[2], orbitcoords[3], orbitcoords[4], orbitcoords[5])
            
            r0 = np.zeros((3,1))
            r0[0, 0], r0[1,0], r0[2,0] = orbitcoords[0], orbitcoords[1], orbitcoords[2]
            y0 = np.array([np.linalg.norm(r0), 0, 0])

            perifocal = rfr.frame_rotation(y0, thetas= [np.degrees(RAAN), np.degrees(inc), np.degrees(w)], euler= 313)

        if coord == 'classical':
            propped_orbit = prop.rk4_2Body(cb= 'earth', dt= dt, tspan = tspan, classical= orbitcoords, show= False)
            t, ys = propped_orbit.propagate_orbit()
            yin = ys[:, :3]

            r0, v0 = cl.class2cart(pd.earth['mu'], orbitcoords[0], orbitcoords[1], orbitcoords[2], orbitcoords[3], orbitcoords[4], orbitcoords[5])
            y0 = np.array([np.linalg.norm(r0), 0, 0])

            perifocal = rfr.frame_rotation(y0, thetas= [np.degrees(orbitcoords[3]), np.degrees(orbitcoords[2]), np.degrees(orbitcoords[4])], euler= 313)

        latlongs, _ , _ = rfr.eci2latlong(yin, np.linspace(0, tspan, steps))

        # Check for ascending Node(s)
        q = 0
        ascends = []
        for step in range(steps): 
            if latlongs[step, 0] >= 0 and latlongs[step-1, 0] < 0:
                ascends.append(q)
            q += 1

        a, e, inc, RAAN, w, M = ct.cart2class(pd.earth['mu'], ys[0, 0], ys[0, 1], ys[0, 2], ys[0, 3], ys[0, 4], ys[0, 5])
        h = np.sqrt(pd.earth['mu'] * a * (1 - e**2))
        p = h**2 / pd.earth['mu']

        # Check for spacecraft location at give time
        q = 0
        if track == True:
            for step in range(steps):
                closest_time = abs(ttt - dt * step)
                if closest_time < dt:
                    sc_step = step
                q += 1
            
        sc_locationX, sc_locationY, sc_locationZ = ys[sc_step, 0], ys[sc_step, 1], ys[sc_step, 2]

        # Storing orbit characteristics to plot against time
        elements = np.zeros((steps, 6))
        velos = np.zeros((steps))
        accels = np.zeros((steps))
        alts = np.zeros((steps))
        fpas = np.zeros((steps))
        nus = np.zeros((steps))
        Es = np.zeros((steps))

        j = 0
        for yn in ys:

            velo = np.linalg.norm(yn[3:])
            velos[j] = velo

            rcurrent = np.linalg.norm(yn[:3])
            ax, ay, az = -(r0 * pd.earth['mu'] / rcurrent**3)
            accel = np.linalg.norm(ax, ay, az)
            accels[j] = accel / 9.8e-3     # Converts from km/s^2 to Gs

            alt = rcurrent - pd.earth['eq']
            alts[j] = alt

            nu = np.arccos((p / rcurrent - 1) / e)
            if np.dot([yn[0], yn[1], yn[2]],[yn[3], yn[4], yn[5]]) < 0: nu = 2*np.pi - nu
            fpa = np.arcsin(pd.earth['mu'] * e * np.sin(nu) / (velo * h))
            fpas[j] = np.degrees(fpa)
            nus[j] = np.degrees(nu)

            E = np.arccos((rcurrent / a - 1) / -e)
            Es[j] = np.degrees(E)

            j += 1

        ##########
        # FRAMES #
        ##########
        # Finding Perifocal and Spacecraft Fixed Frame
        peri = perifocal.find_rotation_matrix()

        nu_current = nus[sc_step]
        if coord == 'cartesian':
            raan = np.degrees(RAAN)
            inclination = np.degrees(inc)
            aop = np.degrees(w)
        else:
            raan = np.degrees(orbitcoords[3])
            inclination = np.degrees(orbitcoords[2])
            aop = np.degrees(orbitcoords[4])

        thet = aop + nu_current
        sc_fixed = rfr.frame_rotation(y0, thetas= [raan, inclination, thet], euler= 313)
        scf = sc_fixed.find_rotation_matrix()

        plt.style.use('dark_background')
        ########################
        # Plotting Earth Orbit #
        ########################
        ax = plt.figure().add_subplot(projection= '3d')

        # Plotting Orbit about Central Body
        ax.plot(ys[:, 0], ys[:, 1], ys[:, 2], color= 'darkmagenta', label= f'{labels[n]} about Earth')
        ax.plot(ys[ascends[0], 0], ys[ascends[0], 1], ys[ascends[0], 2], 'r*', markersize= 12, label= f'Ascending Node(s)')

        if track == True:
            ax.plot(sc_locationX, sc_locationY, sc_locationZ, 'c^', markersize= 12, label= f'Spacecraft Location after {round(ttt * time_conversion, 2)} {time_units}')     # Current position of s/c
            ax.plot(ys[0, 0], ys[0, 1], ys[0, 2], 'y^', markersize= 12, label= f'Initial Position of Spacecraft')     # Initial position of s/c)

        # Plotting Central Body
        _u,_v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
        _x = pd.earth['eq'] * np.cos(_u)*np.sin(_v)
        _y = pd.earth['eq'] * np.sin(_u)*np.sin(_v)
        _z = pd.earth['eq'] * np.cos(_v)

        ax.plot_surface(_x,_y,_z, cmap= pd.earth['clrmap'], alpha= 0.45)

        l = y0[0]
        bigl = l * 1.1
        x, y, z = [[bigl,0,0], [0,bigl,0], [0,0,bigl]]
        zeros = [[0,0,0], [0,0,0], [0,0,0]]

        xprime, yprime, zprime = peri[0, :] * l, peri[1, :] * l, peri[2, :] * l
        rhat, thetahat, hhat = scf[0, :] * l, scf[1, :] * l, scf[2, :] * l
        #sc_zeros = [[sc_locationX, 0, 0], [0, sc_locationY, 0], [0, 0, sc_locationZ]]
        #sc_zeros = [[sc_locationX, sc_locationY, sc_locationZ], [0, 0, 0], [0, 0, 0]]
        adjust = 0.5*l
        sc_zeros = [[sc_locationX, sc_locationX, sc_locationX], [sc_locationY, sc_locationY, sc_locationY], [sc_locationZ, sc_locationZ, sc_locationZ]]
        #sc_zeros = [[0, sc_locationX, sc_locationX], [sc_locationY, 0, sc_locationY], [sc_locationZ, sc_locationZ, 0]]
        #sc_zeros = [[sc_locationX, sc_locationX, sc_locationX], [0,0,0], [0,0,0]]  # Kinda Works ???

        # Plotting Axes #
        ax.quiver(zeros, zeros, zeros, x, y, z , color= 'cornflowerblue', label= 'Earth-Centered Inertial Frame (J2000)')    # Original x-y-z axes
        ax.text( 1.5*l, 0, 0, 'X', color = 'cornflowerblue' )
        ax.text( 0, 1.5*l, 0, 'Y', color = 'cornflowerblue'  )
        ax.text( 0, 0, 1.5*l, 'Z', color = 'cornflowerblue' )

        if disp_peri == True:
            ax.quiver(zeros, zeros, zeros, xprime, yprime, zprime, color= 'lime', label= 'Perifocal Frame')
            ax.text( peri[ 0, 0 ] * 1.2*l, peri[ 1, 0 ] * 1.5*l, peri[ 2, 0 ] * 1.5*l, "X'",
                    color = 'lime' )
            ax.text( peri[ 0, 1 ] * 1.2*l, peri[ 1, 1 ] * 1.5*l, peri[ 2, 1 ] * 1.5*l, "Y'",
                color = 'lime' )
            ax.text( peri[ 0, 2 ] * 1.2*l, peri[ 1, 2 ] * 1.5*l, peri[ 2, 2 ] * 1.5*l, "Z'",
                color = 'lime' )

        if track == True:
            ax.quiver(sc_locationX, sc_locationY, sc_locationZ, 0.5*rhat, 0.5*thetahat, 0.5*hhat, color= 'magenta',  label= 'Spacecraft Body-Fixed Frame')
            ax.text( 1.65*sc_locationX, 1.65*sc_locationY, 1.65*sc_locationZ, "r",
                color = 'magenta', )


        plt.xlim(-l*1.35, l*1.35)
        plt.ylim(-l*1.35, l*1.35)
        ax.set_zlim(-l*1.35, l*1.35)

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.legend(bbox_to_anchor=(0.85, 0.9), loc='lower center')

        #########################
        # Plotting Groundtracks #
        #########################

        plt.figure()
        cs= ['C9', 'c', 'C5', 'C1', 'C6', 'w', 'b', 'C6']

        # Loading Coastline latlongs
        coast_coords = np.genfromtxt('coastlines.csv', delimiter= ',')
        # Plotting Coastlines
        plt.plot(coast_coords[:, 0], coast_coords[:, 1], 'o', markersize= 0.3, color= 'magenta')

        # Plotting Orbits

        plt.plot(latlongs[0, 1], latlongs[0, 0], cs[n] + 'o', label= labels[n])
        plt.plot(latlongs[1:, 1], latlongs[1:, 0], cs[n] + 'o', markersize= 0.75)

        # Plotting Ascending Node(s)   
        count = 0
        for ascend in ascends:
            if count < 1:
                plt.plot(latlongs[ascend, 1], latlongs[ascend, 0], 'r*', markersize= 14, label= 'Ascending Node(s)')
            else:
                plt.plot(latlongs[ascend, 1], latlongs[ascend, 0], 'r*', markersize= 14)
            
            count += 1

        # Plotting S/C Location
        if track == True:
            plt.plot(latlongs[sc_step, 1], latlongs[sc_step, 0], 'c^', markersize= 12, label= f'Spacecraft Location after {round(ttt * time_conversion, 2)} {time_units}')
            plt.plot(latlongs[0, 1], latlongs[0, 0], 'g^', markersize= 12, label= 'Initial Position of Spacecraft')

        # Plotting Cities
        cities = gt.city_dict()
        city_names = gt.city_list1
        i = 0
        for city in city_names:
            citycoords = cities[city]

            citycolor = 'white'
            if i <= 4:
                citycolor = 'lime'

            plt.plot([citycoords[1]], [citycoords[0]], 'o', markersize= 3, color= 'palegreen')

            # Alternating annotation above or below city
            if i % 2 == 0:
                xytext = (0, 3)
            else: 
                xytext = (0, -10)
            
            plt.annotate(city, [citycoords[1], citycoords[0]],
                textcoords= 'offset points', xytext= xytext, ha= 'center', color= citycolor, fontsize= 'small')
            
            i += 1
            
        plt.grid(linestyle= 'dotted')
        plt.xlabel('Longitude (degrees $^\circ$)')
        plt.ylabel('Latitude (degrees $^\circ$)')
        plt.xlim(-180, 180)
        plt.ylim(-90, 90)
        plt.legend(bbox_to_anchor=(0.825, 1), loc='lower center')

        ###################################
        # Plotting Flight Characteristics #
        ###################################

        fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2)
        t = t * time_conversion

        fig.suptitle(f'{labels[n]} Properties')

        ax1.plot(t, velos, color= 'red')
        ax1.set_xlabel(f'Time ({time_units})')
        ax1.set_ylabel('Velocity (km/s)')

        ax2.plot(t, alts, color= 'cyan')
        ax2.set_xlabel(f'Time ({time_units})')
        ax2.set_ylabel('Altitude (km)')

        ax3.plot(t, accels, color= 'magenta')
        ax3.set_xlabel(f'Time ({time_units})')
        ax3.set_ylabel('Acceleration (Gs)')

        ax4.plot(t, nus, color= 'goldenrod')
        ax4.set_xlabel(f'Time ({time_units})')
        ax4.set_ylabel('True Anomaly (degrees $^\circ$)')

        ax5.plot(t, fpas, color= 'lime')
        ax5.set_xlabel(f'Time ({time_units})')
        ax5.set_ylabel('Flight Path Angle (degrees $^\circ$)')

        ax6.plot(t, Es, color= 'blueviolet')
        ax6.set_xlabel(f'Time ({time_units})')
        ax6.set_ylabel('Eccentric Anomaly (degrees $^\circ$)')

        # Printing results
        if track == True:

            root = tk.Tk()     # Creates window of application
            root.geometry('1250x500+225+225')   # Adjusts size of window
            root.title('Orbit Propagation Results')
            mainframe = tk.Frame(root, background='black')    # Places widgets inside a frame (makes more complex GUIs easier to work with)
            mainframe.pack(fill='both', expand='True')

            general_label1 = ttk.Label(mainframe, text= f'For {labels[n]},',
                                    foreground='white', background='black', font=("Brass Mono", 21))
            general_label1.grid(row=0, column=0, padx=10, pady=15)
            general_label2 = ttk.Label(mainframe, text= f'when the spacecraft is',
                                    foreground='white', background='black', font=("Brass Mono", 21))
            general_label2.grid(row=1, column=0, padx=0, pady=15)
            general_label3 = ttk.Label(mainframe, text= f'{round(ttt * time_conversion, 2)} {time_units} into its journey...',
                                    foreground='white', background='black', font=("Brass Mono", 21))
            general_label3.grid(row=1, column=1, padx=0, pady=15)
            # Displaying Velocity
            velo_label1 = ttk.Label(mainframe, text= f'Velocity =', foreground='red', background='black', font=("Brass Mono", 15))
            velo_label1.grid(row=2, column=0, pady=15)
            velo_label2 = ttk.Label(mainframe, text= f'{round(velos[sc_step], 2)} km/s', foreground='white', background='black', font=("Brass Mono", 15))
            velo_label2.grid(row=2, column=1, pady=15)
            # Displaying Altitude
            alt_label1 = ttk.Label(mainframe, text= f'Altitude =', foreground='cyan', background='black', font=("Brass Mono", 15))
            alt_label1.grid(row=3, column=0, pady=15)
            alt_label2 = ttk.Label(mainframe, text= f'{round(alts[sc_step], 2)} km', foreground='white', background='black', font=("Brass Mono", 15))
            alt_label2.grid(row=3, column=1, pady=15)
            # Displaying Acceleration
            accel_label1 = ttk.Label(mainframe, text= f'Acceleration =', foreground='magenta', background='black', font=("Brass Mono", 15))
            accel_label1.grid(row=4, column=0, pady=15)
            accel_label2 = ttk.Label(mainframe, text= f'{round(accels[sc_step]  * 9.8, 2)} km/s^2',
                                     foreground='white', background='black', font=("Brass Mono", 15))
            accel_label2.grid(row=4, column=1, pady=15)
            accel_label3 = ttk.Label(mainframe, text= f'{round(accels[sc_step], 2)} Gs',
                                     foreground='white', background='black', font=("Brass Mono", 15))
            accel_label3.grid(row=4, column=2, pady=15)
            # Displaying True Anomaly
            nu_label1 = ttk.Label(mainframe, text= f'True Anomaly =', foreground='goldenrod', background='black', font=("Brass Mono", 15))
            nu_label1.grid(row=5, column=0, pady=15)
            nu_label2 = ttk.Label(mainframe, text= f'{round(nus[sc_step], 2)} degrees', foreground='white', background='black', font=("Brass Mono", 15))
            nu_label2.grid(row=5, column=1, pady=15)
            # Displaying Flight Path Angle
            fpa_label1 = ttk.Label(mainframe, text= f'Flight Path Angle =', foreground='lime', background='black', font=("Brass Mono", 15))
            fpa_label1.grid(row=6, column=0, pady=15)
            fpa_label2 = ttk.Label(mainframe, text= f'{round(fpas[sc_step], 2)} degrees', foreground='white', background='black', font=("Brass Mono", 15))
            fpa_label2.grid(row=6, column=1, pady=15)
            # Display Eccentric Anomaly
            E_label1 = ttk.Label(mainframe, text= f'Eccentric Anomaly =', foreground='blueviolet', background='black', font=("Brass Mono", 15))
            E_label1.grid(row=7, column=0, pady=15)
            E_label2 = ttk.Label(mainframe, text= f'{round(Es[sc_step], 2)} degrees', foreground='white', background='black', font=("Brass Mono", 15))
            E_label2.grid(row=7, column=1, pady=15)


        n += 1

    plt.show()
