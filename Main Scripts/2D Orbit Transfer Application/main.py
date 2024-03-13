# ---------------------- 2D Circular Orbit Transfer w/ GUI --------------------- #

import tkinter as tk
from tkinter import *
from tkinter import ttk
from PIL import ImageTk, Image

import math as m
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.animation import PillowWriter
import numpy as np
import planorbchar as oc
import choosespacecraft as sp
import time
import warnings

class TransferSim():
    def __init__(self):
        self.root = tk.Tk()     # Creates window of application
        self.root.geometry('1000x1000+225+25')   # Adjusts size of window
        self.root.title('2D Orbit Transfer')
        self.mainframe = tk.Frame(self.root, background='black')    # Places widgets inside a frame (makes more complex GUIs easier to work with)
        self.mainframe.pack(fill='both', expand='True')

        self.orbit_icon = tk.PhotoImage(file= "2d_transfer_Pics/Orbit Icon.png") # Setting Widget Icon
        self.root.iconphoto(False, self.orbit_icon)

        self.unit_label = ttk.Label(self.mainframe, text='Transfer Type: ', foreground='white', background='black', font=("Brass Mono", 15))
        self.unit_label.grid(row=0, column=0, pady=15)
        self.input_label = ttk.Label(self.mainframe, text='Spacecraft: ', foreground='white', background='black', font=("Brass Mono", 15))
        self.input_label.grid(row=1, column=0, pady=15)     
        set_init_param_button = ttk.Button(self.mainframe, text='Set Initial Parameters', command=self.initialize_parameters, width=21)
        set_init_param_button.grid(row=2, column=0, padx=10, pady=10)

        self.depart_planet_label = ttk.Label(self.mainframe, text='Departure Planet: ', foreground='white', background='black', font=("Brass Mono", 15))
        self.depart_planet_label.grid(row=3, column=0, pady=15)
        self.arrive_planet_label = ttk.Label(self.mainframe, text='Arrival Planet: ', foreground='white', background='black', font=("Brass Mono", 15))
        self.arrive_planet_label.grid(row=4, column=0, pady=15)

        self.depart_alt_label = ttk.Label(self.mainframe, text='Departure Altitude (km): ', foreground='white', background='black', font=("Brass Mono", 15))
        self.depart_alt_label.grid(row=3, column=2, padx=5, pady=15)
        self.arrive_alt_label = ttk.Label(self.mainframe, text='Arrival Altitude (km): ', foreground='white', background='black', font=("Brass Mono", 15))
        self.arrive_alt_label.grid(row=4, column=2, padx=5, pady=15)

        set_planets_button = ttk.Button(self.mainframe, text='Set Planets', command=self.set_planets, width=21)
        set_planets_button.grid(row=6, column=0, padx=5, pady=10)

        self.type_label = ttk.Label(self.mainframe, text='Orbit Type: ', foreground='white', background='black', font=("Brass Mono", 15))
        self.type_label.grid(row=5, column=0, pady=10)

        ########################
        # Entry/Dropdown Boxes #
        ########################

        self.transfer_types = ['Single Planet', 'Two Planet', 'Planet to Moon', 'Moon to Planet']
        self.bodies_list = ['Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto']
        self.moon_planet_list = ['Earth-Moon', 'Saturn-Titan']
        self.orbit_types = ['Hohmann', 'Bi-Elliptic']

        self.set_transfer_type = ttk.Combobox(self.mainframe, values=self.transfer_types)
        self.set_transfer_type.grid(row=0, column=1, sticky='NWES', pady=10)
        self.departure_list = []
        self.arrival_list = []

        self.spacecraft_options = ['Voyager', 'Orion', 'Starship', 'Falcon 9']
        self.set_spacecraft = ttk.Combobox(self.mainframe, values=self.spacecraft_options)
        self.set_spacecraft.grid(row=1, column=1, sticky='NWES', pady=10)

        self.set_depart_planet = ttk.Combobox(self.mainframe, values=self.departure_list)
        self.set_depart_planet.grid(row=3, column=1, sticky='NWES', pady=10)
        self.set_arriv_planet = ttk.Combobox(self.mainframe, values=self.arrival_list)
        self.set_arriv_planet.grid(row=4, column=1, sticky='NWES', pady=10)

        self.set_depart_alt = ttk.Entry(self.mainframe)
        self.set_depart_alt.grid(row=3, column=3, sticky='NWES', pady=10)
        self.set_arriv_alt = ttk.Entry(self.mainframe)
        self.set_arriv_alt.grid(row=4, column=3, sticky='NWES', pady=10)

        self.set_orbit_type = ttk.Combobox(self.mainframe, values=self.orbit_types)
        self.set_orbit_type.grid(row=5, column=1, sticky='NWES', pady=10)

        self.set_go_button = ttk.Button(self.mainframe, text='Go', command=self.run_sim, width=21)
        self.set_go_button.grid(row=7, column=0, padx=10, pady=10)

        ##########
        # Images #
        ##########

        # Initial Parameters
        self.trsf_type = 0

        self.root.mainloop()    # Dself.isplays the window on screen

    def initialize_parameters(self):

        self.rad = np.pi/180  #self.radian conversion from degrees
        self.g = 0.00981  #Acceleration due to gravity on Earth in km/s^2

        #Characteristics of the Sun
        [self.muS,eqS,MS,self.clrS] = oc.planchar('Sun')
        #Characteristics of Earth's Orbit around the Sun
        [a_eth,e_eth,rp_eth,ra_eth,b_eth,p_eth,h_eth] = oc.orbchar('Earth')

        sc = self.set_spacecraft.get()
        [self.sc_type,self.scclr,self.oclr,self.m_sc,self.isp] = sp.scchsr(sc)

        if self.set_transfer_type.get() == 'Single Planet':
            self.trsf_type = '1'
            self.set_depart_planet.config(values=self.bodies_list)
            self.set_arriv_planet.config(values='N/A')
        elif self.set_transfer_type.get() == 'Two Planet':
            self.trsf_type = '2'
            self.set_depart_planet.config(values=self.bodies_list)
            self.set_arriv_planet.config(values=self.bodies_list)
        elif self.set_transfer_type.get() == 'Planet to Moon':
            self.trsf_type = 'pm'
            self.set_depart_planet.config(values=self.moon_planet_list)
            self.set_arriv_planet.config(values='N/A')
        elif self.set_transfer_type.get() == 'Moon to Planet':
            self.trsf_type = 'mp'
            self.set_depart_planet.config(values=self.moon_planet_list)
            self.set_arriv_planet.config(values='N/A')



    def set_planets(self):

        ###################################################
        # Deleting Residual Labels for Bi-Elliptic Inputs #
        ###################################################

        try:
            self.int_rad_label.destroy()
            self.set_int_rad.destroy()
            self.int_rad_exp_label.destroy()
            self.set_int_rad_exp.destroy()
            self.int_rad_km_label.destroy()
            
        except:
            pass

        self.set_go_button.grid(row=7, column=0, padx=10, pady=10)

        ############################
        # Basic Orbit Calculations #
        ############################

        if self.trsf_type == '1':

            #For transfer between two orbits in a SINGLE CELESTIAL BODY'S Sphere of Influence
            self.planet = self.set_depart_planet.get()
            #Determining this body's characteristics
            [self.mu,eq,M,self.clr] = oc.planchar(self.planet)
            self.mu1 = self.mu
            self.mu2 = self.mu

            if self.planet == 'Earth':

                self.sc_start = self.set_depart_alt.get()
            
                if self.sc_start == 'LEO' or self.sc_start == 'leo':
                    self.alt1 = 160
                    self.r1 = eq + self.alt1
                
                elif self.sc_start == 'GEO' or self.sc_start == 'geo':
                    self.alt1 = 35789
                    self.r1 = eq + self.alt1

                else:
                    #Determining self.r1
                    self.alt1 = float(self.set_depart_alt.get())
                    self.r1 = eq + self.alt1
            
            else:
                #Determining self.r1
                self.alt1 = float(self.set_depart_alt.get())
                self.r1 = eq + self.alt1

            #Characteristics of initial orbit (Circular)
            e1 = 0 #Because the spacecraft's orbit is assumed circular 
            a1 = self.r1/(1+e1)
            b1 = a1
            p1 = a1*(1-e1**2)
            self.h1 = m.sqrt(p1*self.mu)

            #Determining self.r2
            self.alt2st = self.set_arriv_alt.get()
            if self.alt2st == 'GEO' or self.alt2st == 'geo':
                self.alt2 = 35786
            elif self.alt2st == 'LEO' or self.alt2st == 'leo':
                    self.alt2 = 160
                    self.r2 = eq + self.alt2
            else:
                self.alt2 = float(self.alt2st)
            self.r2 = eq + self.alt2
            #Characteristics of final orbit (Circular)
            e2 = 0 #Because the spacecraft's orbit is assumed circular 
            a2 = self.r2/(1+e2)
            b2 = a2
            p2 = a2*(1-e2**2)
            self.h2 = m.sqrt(p2*self.mu)
            
        #For transfer between TWO DIFFERENT BODIES' ORBITS
        elif self.trsf_type == '2':

            self.mu = self.muS  #Used when referencing gravitational parameter on transfer orbit
            
            #For transfer between TWO DIFFERENT CELESTIAL BODIES
            self.planet1 = self.set_depart_planet.get()
            self.planet2 = self.set_arriv_planet.get()
            #Determining these bodies' characteristics
            [self.mu1,eq1,M1,self.clr1] = oc.planchar(self.planet1)
            [self.mu2,eq2,M2,self.clr2] = oc.planchar(self.planet2)
            
            #Characteristics of body 1's orbit around the Sun
            [a_b1,e_b1,rp_b1,ra_b1,b_b1,p_b1,h_b1] = oc.orbchar(self.planet1)
            #Circularizing Orbit 1 and determining self.r1
            rb1 = a_b1  #Because e=0 for a circular orbit
            #self.alt1 = float(input(f'Enter the altitude above {self.planet1} (in km) from which \nthe spacecraft will depart: '))
            self.alt1 = float(self.set_depart_alt.get())
            self.r1 = rb1 + eq1 + self.alt1

            #Characteristics of body 2's orbit around the Sun
            [a_b2,e_b2,rp_b2,ra_b2,b_b2,p_b2,h_b2] = oc.orbchar(self.planet2)
            #Circularizing Orbit 2 and determining self.r2
            rb2 = a_b2  #Because e=0 for a circular orbit
            #self.alt2 = float(input(f'Enter the altitude above {self.planet2} (in km) from which \nthe spacecraft will arrive: '))
            self.alt2 = float(self.set_arriv_alt.get())
            self.r2 = rb2 + eq2 + self.alt2

            #Characteristics of initial orbit (Circular)
            e1 = 0  #Because the spacecraft's orbit is assumed circular 
            a1 = self.r1/(1+e1)
            b1 = a1
            p1 = a1*(1-e1**2)
            self.h1 = m.sqrt(p1*self.mu1)
            #Characteristics of final orbit (Circular)
            e2 = 0  #Because the spacecraft's orbit is assumed circular 
            a2 = self.r2/(1+e2)
            b2 = a2
            p2 = a2*(1-e2**2)
            self.h2 = m.sqrt(p2*self.mu2)

        #elif self.trsf_type == 'pm' or self.trsf_type == 'mp':
        else:
            
            moonset = self.set_depart_planet.get()

            if moonset == 'Earth-Moon':
                whichmoon = 'E'
            elif moonset == 'Saturn-Titan':
                whichmoon = 'T'


            if (whichmoon == 'E' or whichmoon == 'e') and (self.trsf_type == 'pm'):
                #Determining Earth and Moon characteristics
                [self.mu1,eq1,M1,self.clr1] = oc.planchar('Earth')
                [self.mu2,eq2,M2,self.clr2] = oc.planchar('Moon')
                [a_mn,e_mn,rp_mn,ra_mn,b_mn,p_mn,h_mn] = oc.orbchar('Moon')
                self.mu = self.mu1
                r_mn = a_mn
                self.body1 = 'Earth'
                self.body2 = 'Moon'

            elif (whichmoon == 'E' or whichmoon == 'e') and (self.trsf_type == 'mp'):
                #Determining Earth and Moon characteristics
                [self.mu2,eq2,M2,self.clr2] = oc.planchar('Earth')
                [self.mu1,eq1,M1,self.clr1] = oc.planchar('Moon')
                [a_mn,e_mn,rp_mn,ra_mn,b_mn,p_mn,h_mn] = oc.orbchar('Moon')
                self.mu = self.mu2
                r_mn = a_mn
                self.body1 = 'Moon'
                self.body2 = 'Earth'

            elif (whichmoon == 'T' or whichmoon == 't') and (self.trsf_type == 'pm'):
                #Determining Earth and Moon characteristics
                [self.mu1,eq1,M1,self.clr1] = oc.planchar('Saturn')
                [self.mu2,eq2,M2,self.clr2] = oc.planchar('Titan')
                [a_mn,e_mn,rp_mn,ra_mn,b_mn,p_mn,h_mn] = oc.orbchar('Titan')
                self.mu = self.mu1
                r_mn = a_mn
                self.body1 = 'Saturn'
                self.body2 = 'Titan'

            elif (whichmoon == 'T' or whichmoon == 't') and (self.trsf_type == 'mp'):
                #Determining Earth and Moon characteristics
                [self.mu1,eq1,M1,self.clr1] = oc.planchar('Titan')
                [self.mu2,eq2,M2,self.clr2] = oc.planchar('Saturn')
                [a_mn,e_mn,rp_mn,ra_mn,b_mn,p_mn,h_mn] = oc.orbchar('Titan')
                self.mu = self.mu2
                r_mn = a_mn
                self.body1 = 'Titan'
                self.body2 = 'Saturn'

            #Determining self.alt1
            self.alt1 = float(self.set_depart_alt.get())        
            #Determing self.alt2
            self.alt2 = float(self.set_arriv_alt.get())

            #Determining self.r1 and self.r2
            if self.trsf_type == 'pm':
                self.r1 = eq1 + self.alt1
                self.r2 = eq2 + a_mn + self.alt2
            elif self.trsf_type == 'mp':
                self.r1 = eq1 + a_mn + self.alt1
                self.r2 = eq2 + self.alt2
            
            #Characteristics of initial orbit (Circular)
            e1 = 0  #Because the spacecraft's orbit is assumed circular 
            a1 = self.r1/(1+e1)
            b1 = a1
            p1 = a1*(1-e1**2)
            self.h1 = m.sqrt(p1*self.mu1)
            #Characteristics of final orbit (Circular)
            e2 = 0  #Because the spacecraft's orbit is assumed circular 
            a2 = self.r2/(1+e2)
            b2 = a2
            p2 = a2*(1-e2**2)
            self.h2 = m.sqrt(p2*self.mu2)

        orb_type = self.set_orbit_type.get()
        if orb_type == 'Bi-Elliptic':

            nu1 = 270*self.rad
            nuI = 90*self.rad
            nu2 = 270*self.rad

            self.r1_rep = self.r1  #Reported semi-major axis
            self.exp_r1 = 0  #Reported exponent for semi-major axis measurement
            for n in range(0,25):
                if not self.r1_rep < 10:
                    self.r1_rep = round(self.r1/(10**n),3)
                    self.exp_r1 = n

            self.int_rad_label = ttk.Label(self.mainframe, text=f'(Set Intermediate Radius\ngreater than {self.r1_rep}^{self.exp_r1} km): ', foreground='white', background='black', font=("Brass Mono", 15))
            self.int_rad_label.grid(row=7, column=0, pady=1)
            self.set_int_rad = ttk.Entry(self.mainframe)
            self.set_int_rad.grid(row=7, column=1, sticky='NWES', pady=1)
            self.int_rad_exp_label = ttk.Label(self.mainframe, text='^', foreground='white', background='black', font=("Brass Mono", 15))
            self.int_rad_exp_label.grid(row=7, column=2, pady=1)
            self.set_int_rad_exp = ttk.Entry(self.mainframe)
            self.set_int_rad_exp.grid(row=7, column=3, sticky='NWES', pady=1)
            self.int_rad_km_label = ttk.Label(self.mainframe, text='km', foreground='white', background='black', font=("Brass Mono", 15))
            self.int_rad_km_label.grid(row=7, column=4, padx=1, pady=1)

            self.set_go_button.grid(row=8, column=0, padx=10, pady=10)      # Moving the go button



    def run_sim(self):

        orb_type = self.set_orbit_type.get()

        ###########################################
        # Planet Images/Non-Numerical Interfacing #
        ###########################################

        if self.trsf_type == '1':

            try:
                self.type_label.destroy()
                self.img1_label.destroy()
                self.img2_label.destroy()
                self.to_label.destroy()
            except:
                pass

            self.type_label = Label(self.mainframe, text=f'Intra-Planetary\n{self.set_orbit_type.get()} Transfer', foreground='White', background='Black', font=("Brass Mono", 17))
            self.type_label.grid(row=9, column=1, sticky='NWES')

            self.img = Image.open(f"2d_transfer_Pics/{self.planet}_Pic.png")
            self.img_resized = self.img.resize((100,100), Image.ANTIALIAS)
            self.new_img = ImageTk.PhotoImage(self.img_resized)
            self.img_label = Label(self.mainframe, image=self.new_img, background='Black')
            self.img_label.grid(row=9, column=2, sticky='NWES')

        else:

            if self.trsf_type != '2':
                planet1 = self.body1
                planet2 = self.body2
            else:
                planet1 = self.planet1
                planet2 = self.planet2

            try:
                self.type_label.destroy()
                self.img_label.destroy()
                self.sc_img_label.destroy()
            except:
                pass

            try:
                self.type_label.destroy()
                self.img1_label.destroy()
                self.img2_label.destroy()
                self.to_label.destroy()
            except:
                pass

            if planet1 == 'Moon' or planet2 == 'Moon':
                self.type_label = Label(self.mainframe, text=f'Lunar\n{self.set_orbit_type.get()} Transfer', foreground='White', background='Black', font=("Brass Mono", 17))
            else:
                self.type_label = Label(self.mainframe, text=f'Interplanetary\n{self.set_orbit_type.get()} Transfer', foreground='White', background='Black', font=("Brass Mono", 17))
            self.type_label.grid(row=9, column=0, sticky='NWES')

            self.img1 = Image.open(f"2d_transfer_Pics/{planet1}_Pic.png")
            self.img1_resized = self.img1.resize((100,100), Image.ANTIALIAS)
            self.new_img1 = ImageTk.PhotoImage(self.img1_resized)
            self.img1_label = Label(self.mainframe, image=self.new_img1, background='Black')
            self.img1_label.grid(row=9, column=1, sticky='NWES')

            self.to_label = Label(self.mainframe, text='→', foreground='White', background='Black', font=("Brass Mono", 30))
            self.to_label.grid(row=9, column=2, sticky='NWES')

            self.img2 = Image.open(f"2d_transfer_Pics/{planet2}_Pic.png")
            self.img2_resized = self.img2.resize((100,100), Image.ANTIALIAS)
            self.new_img2 = ImageTk.PhotoImage(self.img2_resized)
            self.img2_label = Label(self.mainframe, image=self.new_img2, background='Black')
            self.img2_label.grid(row=9, column=3, sticky='NWES')

        if orb_type == 'Hohmann' or orb_type == 'hm' or orb_type == 'Hm':

            ############################
            # Deleting Residual Labels #
            ############################

            self.set_go_button.grid(row=7, column=0, padx=10, pady=10)

            try:

                self.print_ellip_opt_1.destroy()
                self.print_ellip_opt_2.destroy()

                self.print_at1_1.destroy()
                self.print_at1_2.destroy()

                self.print_et1_1.destroy()
                self.print_et1_2.destroy()

                self.print_at2_1.destroy()
                self.print_at2_2.destroy()

                self.print_et2_1.destroy()
                self.print_et2_2.destroy()
                
                self.print_mpT_text1.destroy()
                self.print_mpT_text2.destroy()

                self.print_MF_text1.destroy()
                self.print_MF_text2.destroy()

                self.sc_name.destroy()
        
                self.print_deltv1_1.destroy()
                self.print_deltv1_2.destroy()

                self.print_deltv2_1.destroy()
                self.print_deltv2_2.destroy()

                self.print_deltv3_1.destroy()
                self.print_deltv3_2.destroy()
            
                self.print_TOF_1.destroy()
                self.print_TOF_2.destroy()
            
            except:
                pass

            try:

                self.print_hohm_opt_1.destroy()
                self.print_hohm_opt_2.destroy()

                self.print_a_1.destroy()
                self.print_a_2.destroy()

                self.print_e_1.destroy()
                self.print_e_2.destroy()

                self.print_mpT_text1.destroy()
                self.print_mpT_text2.destroy()

                self.print_MF_text1.destroy()
                self.print_MF_text2.destroy()

                self.sc_name.destroy()

                self.print_deltv1_1.destroy()
                self.print_deltv1_2.destroy()

                self.print_deltv2_1.destroy()
                self.print_deltv2_2.destroy()

                self.print_TOF_1.destroy()
                self.print_TOF_2.destroy()

            except:
                pass


            ##########################################
            # Defining Parameters for Transfer Orbit #
            ##########################################
            
            nu1 = 270*self.rad
            nu2 = 90*self.rad

            if self.r1 < self.r2:
                rp = self.r1
                ra = self.r2
            elif self.r1 > self.r2:
                rp = self.r2
                ra = self.r1

            #Characteristics of the transfer orbit
            a = 1/2*(ra+rp)
            e = 1-rp/a
            b = a*m.sqrt(1-e**2)
            p = a*(1-e**2)
            h = m.sqrt(self.mu*p)

            ###Finding instantaneous velocity changes###

            #S/C velocities on initial orbit
            self.vr1 = 0  #Spacecraft velocity is entirely in theta hat direction
            vtet1 = self.h1/self.r1 
            v1 = m.sqrt(self.mu1/self.r1) #Velocity magnitude anywhere on LEO circular orbit

            #S/C velocities at perigee on transfer orbit
            #nup = m.acos(1/e*(p/rp-1))
            nup = 0
            nup_deg = nup*180/m.pi
            vrp = self.mu*e*m.sin(nup)/h
            vtetp = h/rp
                
            #Delta V1
            v1_vecr = vrp - self.vr1  #Impulsive burn from LEO to transfer orbit
            v1_vectet = vtetp - vtet1
            v1 = m.sqrt(v1_vecr**2+v1_vectet**2)  #Magnitude of Delta V1

            #S/C velocities on final orbit
            self.vr2 = 0  #Spacecraft velocity is entirely in theta hat direction
            vtet2 = self.h2/self.r2 
            v2 = m.sqrt(self.mu2/self.r2) #Velocity magnitude anywhere on the final circular orbit

            #S/C velocities at apogee on transfer orbit
            #nua = m.acos(1/e*(p/ra-1))
            nua = m.pi
            nua_deg = nua*180/m.pi
            vra = self.mu*e*m.sin(nua)/h
            vteta = h/ra

            #Delta V2
            v2_vecr = self.vr2 - vra  #Impulsive burn from transfer orbit to final orbit
            v2_vectet = vtet2 - vteta
            v2 = m.sqrt(v2_vecr**2+v2_vectet**2)  #Magnitude of Delta V2

            #Calculating Time of Flight in MONTHS AND DAYS (IN USE)
            TOFmd = m.pi*m.sqrt(a**3/self.mu)/3600/24/30 #Approximately calculates time of flight in months
            TOFm = int(TOFmd//1) #Calculates TOF to nearest month
            TOFr1 = TOFmd%1  #Calculates remainder of above TOF calculation
            TOFd = int(TOFr1*30) #Finds number of days left in TOF calculation
                
            #Calculating Time of Flight in HOURS AND MINUTES
            TOFhm = m.pi*m.sqrt(a**3/self.mu)/3600  #Exactly calculates time of flight in hours
            TOFh = int(TOFhm//1)  #Calculates TOF to nearest hour
            TOFr2 = TOFhm%1 #Calculates remainder of above TOF calculation
            TOFmn = int(TOFr2*60)  #Finds number of minutes left in TOF calculation

            if TOFm == 0 and TOFd < 4:
                TOF1 = TOFh
                TOF2 = TOFmn
                unit1 = 'hours'
                unit2 = 'minutes'
            else:
                TOF1 = TOFm
                TOF2 = TOFd
                unit1 = 'months'
                unit2 = 'days'

            #Calculating mass of propellant needed for mission (all in kg)
            mf2 = self.m_sc  #The spacecraft's final mass after burn 2 (equal to dry mass)
            mp2 = mf2*(m.e**(v2/(self.isp*self.g))-1) #Propellant mass needed to execute burn 2
            mf1 = mp2 + mf2  #Propellant mass needed after burn 1
            mp1 = mf1*(m.e**(v1/(self.isp*self.g))-1) #Propellant mass need to execute burn 1
            mpT = mp1 + mp2  #Total propellant mass need for whole mission
            mo = mpT + self.m_sc  #Total mass of S/C at beginning of mission (throw mass)
            MF = mf2/mo*100  #Mass Fraction as a percentage (spacecraft's final mass compared to spacecraft's initial mass)

            #Printing results for Hohmann Transfer

            hohm_opt = ''

            if (1 < self.r2/self.r1 and self.r2/self.r1 < 11.94 or 1 < self.r1/self.r2 and self.r1/self.r2 < 11.94) and MF > 10:
                hohm_opt_1 = "The Hohmann ellipse" 
                hohm_opt_2 = "IS optimal for this transfer."
            else:
                hohm_opt_1 = "The Hohmann ellipse is" 
                hohm_opt_2 = "NOT optimal for this transfer."

            self.print_hohm_opt_1 = ttk.Label(self.mainframe, text=hohm_opt_1, foreground='white', background='black', font=("Brass Mono", 15))
            self.print_hohm_opt_1.grid(row=10, column=0, pady=15)
            self.print_hohm_opt_2 = ttk.Label(self.mainframe, text=hohm_opt_2, foreground='white', background='black', font=("Brass Mono", 15))
            self.print_hohm_opt_2.grid(row=10, column=1, pady=15)

            #Using a more reader-friendly value for semi-major axis and mass of propellant
            a_rep = a  #Reported semi-major axis
            mpT_rep = mpT  #Reported mass of propellant required
            exp_a = 0  #Reported exponent for semi-major axis measurement
            exp_mp = 0  #Reported exponent for mass propellant measurement
            for n in range(0,25):
                if not a_rep < 10:
                    a_rep = a/(10**n)
                    exp_a = n
                if not mpT_rep < 10:
                    mpT_rep = mpT/(10**n)
                    exp_mp = n
                    
            a_text1 = "Semi-Major Axis (a):"
            a_text2 = f"{round(a_rep,3)}*10^{exp_a} km"
            self.print_a_1 = ttk.Label(self.mainframe, text=a_text1, foreground='white', background='black', font=("Brass Mono", 14))
            self.print_a_1.grid(row=11, column=0, pady=1)
            self.print_a_2 = ttk.Label(self.mainframe, text=a_text2, foreground='white', background='black', font=("Brass Mono", 15))
            self.print_a_2.grid(row=12, column=0, pady=1)

            e_text1 = "Eccentricity (e):"
            e_text2 = f"{round(e,3)}"
            self.print_e_1 = ttk.Label(self.mainframe, text=e_text1, foreground='white', background='black', font=("Brass Mono", 14))
            self.print_e_1.grid(row=11, column=1, pady=1)
            self.print_e_2 = ttk.Label(self.mainframe, text=e_text2, foreground='white', background='black', font=("Brass Mono", 15))
            self.print_e_2.grid(row=12, column=1, pady=1)
            
            deltv1_text1 = "ΔV1:"
            deltv1_text2 = f"{round(v1,2)} km/s"
            self.print_deltv1_1 = ttk.Label(self.mainframe, text=deltv1_text1, foreground='white', background='black', font=("Brass Mono", 15))
            self.print_deltv1_1.grid(row=11, column=2, pady=1)
            self.print_deltv1_2 = ttk.Label(self.mainframe, text=deltv1_text2, foreground='white', background='black', font=("Brass Mono", 15))
            self.print_deltv1_2.grid(row=12, column=2, pady=1)

            deltv2_text1 = "ΔV2:"
            deltv2_text2 = f"{round(v2,2)} km/s"
            self.print_deltv2_1 = ttk.Label(self.mainframe, text=deltv2_text1, foreground='white', background='black', font=("Brass Mono", 15))
            self.print_deltv2_1.grid(row=11, column=3, pady=1)
            self.print_deltv2_2 = ttk.Label(self.mainframe, text=deltv2_text2, foreground='white', background='black', font=("Brass Mono", 15))
            self.print_deltv2_2.grid(row=12, column=3, pady=1)

            TOF_text1 = "Time of Flight (TOF):"
            TOF_text2 = f"{round(TOF1,2)} {unit1} and {round(TOF2,2)} {unit2}"
            TOF_text3 = f"and {round(TOF2,2)} {unit2}"
            self.print_TOF_1 = ttk.Label(self.mainframe, text=TOF_text1, foreground='white', background='black', font=("Brass Mono", 14))
            self.print_TOF_1.grid(row=13, column=1, pady=30)
            self.print_TOF_2 = ttk.Label(self.mainframe, text=TOF_text2, foreground='white', background='black', font=("Brass Mono", 15))
            self.print_TOF_2.grid(row=13, column=2, pady=30)

            self.sc_label = Label(self.mainframe, text=f'Spacecraft Used:', foreground='White', background='Black', font=("Brass Mono", 16))
            self.sc_label.grid(row=15,column=0,sticky='NWES')
            
            self.sc_img = Image.open(f"2d_transfer_Pics/{self.set_spacecraft.get()}_Pic.png")
            self.sc_img_resized = self.sc_img.resize((100,100), Image.ANTIALIAS)
            self.new_sc_img = ImageTk.PhotoImage(self.sc_img_resized)
            self.sc_img_label = Label(self.mainframe, image=self.new_sc_img, background='Black')
            self.sc_img_label.grid(row=15, column=1, sticky='NWES')

            self.sc_name = Label(self.mainframe, text=f'{self.set_spacecraft.get()}', foreground='White', background='Black', font=("Brass Mono", 15))
            self.sc_name.grid(row=16,column=1,sticky='NWES')

            mpT_text1 = "Propellant Required:"
            mpT_text2 = f"{round(mpT_rep,3)}*10^{exp_mp} kg"
            self.print_mpT_text1 = ttk.Label(self.mainframe, text=mpT_text1, foreground='white', background='black', font=("Brass Mono", 14))
            self.print_mpT_text1.grid(row=15, column=2, pady=1)
            self.print_mpT_text2 = ttk.Label(self.mainframe, text=mpT_text2, foreground='white', background='black', font=("Brass Mono", 15))
            self.print_mpT_text2.grid(row=16, column=2, pady=1)

            MF_text1 = "Propellant Mass Fraction (MF):"
            MF_text2 = f"{round(MF,2)}%"
            self.print_MF_text1 = ttk.Label(self.mainframe, text=MF_text1, foreground='white', background='black', font=("Brass Mono", 14))
            self.print_MF_text1.grid(row=15, column=3, pady=1)
            self.print_MF_text2 = ttk.Label(self.mainframe, text=MF_text2, foreground='white', background='black', font=("Brass Mono", 15))
            self.print_MF_text2.grid(row=16, column=3, pady=1)

            #Plotting the spacecraft's elliptical orbit
            t = np.linspace(-2*m.pi,2*m.pi,10**6)
            z = np.linspace(0,m.pi,10**6)
            zan = np.linspace(0,np.pi,100)  #Sets values for animation functionto iterate thru
            zz = np.linspace(m.pi,2*m.pi,10**6)

            if self.r1 < self.r2:
                x = (a*e+rp)*np.cos(z)-1/2*(ra-rp)
                xan = (a*e+rp)*np.cos(zan)-1/2*(ra-rp)  #Sets values for animation functionto iterate thru
                xx = (a*e+rp)*np.cos(zz)-1/2*(ra-rp)
            elif self.r1 > self.r2:
                x = (a*e+rp)*np.cos(z)+1/2*(ra-rp)
                xan = (a*e+rp)*np.cos(zan)+1/2*(ra-rp)  #Sets values for animation functionto iterate thru
                xx = (a*e+rp)*np.cos(zz)+1/2*(ra-rp)

            y = b*np.sin(z) 
            yan = b*np.sin(zan)  #Sets values for animation functionto iterate thru
            yy = b*np.sin(zz)
            mark = 8  

            plt.style.use('dark_background')
            fig, ax = plt.subplots(figsize=(9,9))
            ax = plt.axes()
            plt.xlabel("Distance (km)")
            plt.ylabel("Distance (km)")
            
            if self.trsf_type == '1':
                plt.title(f"{self.sc_type} Spacecraft's Hohmann Transfer \n\nfrom {self.planet} Orbit at Altitude of {self.alt1} km to Altitude of {self.alt2} km")
                #Central Body
                ax.plot(0,0,marker='o',
                    markersize = 25, color = self.clr)

            elif self.trsf_type == '2':
                plt.title(f"{self.sc_type} Spacecraft's Hohmann Transfer \nfrom {self.planet1} to {self.planet2}")
                #Sun
                ax.plot(0,0, marker='o',
                    markersize = 25, color='darkorange')
            elif self.trsf_type == 'pm' or self.trsf_type == 'mp':
                plt.title(f"{self.sc_type} Spacecraft's Hohmann Transfer \nfrom {self.body1} to {self.body2}")
                #Planet
                if self.trsf_type == 'pm':
                    color = self.clr1
                elif self.trsf_type == 'mp':
                    color = self.clr2
                ax.plot(0,0, marker='o',
                    markersize = 25, color=color)

            #Transfer Orbit
            ax.plot(x,y,label = "Transfer Orbit",color = 'red')
            ax.plot(xx,yy,linestyle='--',linewidth = 0.4,color = 'red')
            #Orbit of Planet 1
            x1 = self.r1*np.cos(t)
            y1 = self.r1*np.sin(t)
            ax.plot(x1,y1, color='lightsteelblue')
            #Orbit of Planet 2
            x2 = self.r2*np.cos(t)
            y2 = self.r2*np.sin(t)
            ax.plot(x2,y2, color='lightsteelblue')
            
            #Spacecraft Locations 1 and 2
            if self.trsf_type == '1':
                ax.plot(-self.r1*np.sin(nu1),0, marker=(4,0,0),
                    markersize = 12, color='yellow',
                        label = 'Initial S/C Position', alpha=0.4)
                ax.plot(-self.r2*np.sin(nu2),0, marker=(4,0,0),
                    markersize = 12, color='lime',
                        label = 'Final S/C Position', alpha=0.4)
            elif self.trsf_type == '2':  #This plots the location of each body on their respective orbtis
                ax.plot(-self.r1*np.sin(nu1),0, marker='o',
                    markersize = 6.378*3, color=self.clr1, label = f'{self.planet1}')
                ax.plot(-self.r2*np.sin(nu2),0, marker='o',
                    markersize = 6.378*3, color=self.clr2, label = f'{self.planet2}')
            elif self.trsf_type == 'pm':
                ax.plot(-self.r1*np.sin(nu1),0, marker=(4,0,0),
                    markersize = 12, color='yellow',
                        label = 'Initial S/C Position', alpha=0.4)
                ax.plot(-self.r2*np.sin(nu2),0, marker='o',
                    markersize = 6.378*3, color=self.clr2, label = self.body2)
            elif self.trsf_type == 'mp':
                ax.plot(-self.r1*np.sin(nu1),0, marker='o',
                    markersize = 6.378*3, color=self.clr1,
                        label = self.body1)
                ax.plot(-self.r2*np.sin(nu2),0, marker=(4,0,0),
                    markersize = 12, color='lime',
                        label = 'Final S/C Position', alpha=0.4)

            #Plotting spacecraft's movement along transfer orbit
            
            trsf, = ax.plot([],[],marker=(4,0,0),
                            color=self.scclr,markersize=12,label=f'{self.sc_type} Spacecraft')
            def animate(i):
                trsf.set_data(xan[i], yan[i])
                return trsf
            anim = FuncAnimation(fig, animate, frames=100, interval=200, repeat=True)
            #anim.save("EarthTransfer_CustomOrbits.gif")
            
            plt.grid(color = 'lightgray', linestyle='--',linewidth = 0.3)
            plt.legend(bbox_to_anchor=(0.95, 0.85), fontsize='16.5', loc='lower center')
            plt.show()


            
        elif orb_type == 'Bi-Elliptic' or orb_type == 'E' or orb_type == 'e':

            #########################################
            # Destroying Potential Residual Widgets #
            #########################################
            
            try:

                self.print_hohm_opt_1.destroy()
                self.print_hohm_opt_2.destroy()

                self.print_a_1.destroy()
                self.print_a_2.destroy()

                self.print_e_1.destroy()
                self.print_e_2.destroy()

                self.print_mpT_text1.destroy()
                self.print_mpT_text2.destroy()

                self.print_MF_text1.destroy()
                self.print_MF_text2.destroy()

                self.sc_name.destroy()

                self.print_deltv1_1.destroy()
                self.print_deltv1_2.destroy()

                self.print_deltv2_1.destroy()
                self.print_deltv2_2.destroy()

                self.print_TOF_1.destroy()
                self.print_TOF_2.destroy()

            except:
                pass

            try:

                self.print_ellip_opt_1.destroy()
                self.print_ellip_opt_2.destroy()

                self.print_at1_1.destroy()
                self.print_at1_2.destroy()

                self.print_et1_1.destroy()
                self.print_et1_2.destroy()

                self.print_at2_1.destroy()
                self.print_at2_2.destroy()

                self.print_et2_1.destroy()
                self.print_et2_2.destroy()
                
                self.print_mpT_text1.destroy()
                self.print_mpT_text2.destroy()

                self.print_MF_text1.destroy()
                self.print_MF_text2.destroy()

                self.sc_name.destroy()
        
                self.print_deltv1_1.destroy()
                self.print_deltv1_2.destroy()

                self.print_deltv2_1.destroy()
                self.print_deltv2_2.destroy()

                self.print_deltv3_1.destroy()
                self.print_deltv3_2.destroy()
            
                self.print_TOF_1.destroy()
                self.print_TOF_2.destroy()
            
            except:
                pass


            ##########################################
            # Defining Parameters for Transfer Orbits #
            ##########################################

            nu1 = 270*self.rad
            #nuI = 90*self.rad
            nu2 = 270*self.rad

            self.r1_rep = self.r1  #Reported semi-major axis
            self.exp_r1 = 0  #Reported exponent for semi-major axis measurement
            for n in range(0,25):
                if not self.r1_rep < 10:
                    self.r1_rep = round(self.r1/(10**n),3)
                    self.exp_r1 = n
            rI_base = self.set_int_rad.get()
            rI_exp = self.set_int_rad_exp.get()
            rI = float(rI_base)*10**float(rI_exp)

            #Characteristics of the FIRST transfer orbit
            rpt1 = self.r1
            rat1 = rI
            at1 = 1/2*(rat1+rpt1)
            et1 = 1-rpt1/at1
            bt1 = at1*m.sqrt(1-et1**2)
            pt1 = at1*(1-et1**2)
            ht1 = m.sqrt(self.mu*pt1)
            
            ###Finding instantaneous velocity changes###

            #S/C velocities on initial orbit
            self.vr1 = 0  #Spacecraft velocity is entirely in theta hat direction
            vtet1 = self.h1/self.r1 
            v1 = m.sqrt(self.mu1/self.r1) #Velocity magnitude anywhere on LEO circular orbit

            #S/C velocities at perigee on transfer orbit 1
            #nup = m.acos(1/e*(p/rp-1))
            nupt1 = 0
            nupt1_deg = nupt1*180/m.pi
            vrpt1 = self.mu*et1*m.sin(nupt1)/ht1
            vtetpt1 = ht1/rpt1
                
            #Delta V1
            v1_vecr = vrpt1 - self.vr1  #Impulsive burn from LEO to transfer orbit
            v1_vectet = vtetpt1 - vtet1
            v1 = m.sqrt(v1_vecr**2+v1_vectet**2)  #Magnitude of Delta V1

            #S/C velocities at apogee on transfer orbit 1
            #nua = m.acos(1/e*(p/ra-1))
            nuat1 = m.pi
            nuat1_deg = nuat1*180/m.pi
            vrat1 = self.mu*et1*m.sin(nuat1)/ht1
            vtetat1 = ht1/rat1

            #Characteristics of the SECOND transfer orbit
            rpt2 = self.r2
            rat2 = rI
            at2 = 1/2*(rat2+rpt2)
            et2 = 1-rpt2/at2
            bt2 = at2*m.sqrt(1-et2**2)
            pt2 = at2*(1-et2**2)
            ht2 = m.sqrt(self.mu*pt2)

            #S/C velocities at apogee on transfer orbit 2
            #nua = m.acos(1/e*(p/ra-1))
            nuat2 = m.pi
            nuat2_deg = nuat2*180/m.pi
            vrat2 = self.mu*et2*m.sin(nuat2)/ht2
            vtetat2 = ht2/rat2

            #Delta V2
            v2_vecr = vrat2 - vrat1  #Impulsive burn from transfer orbit to final orbit
            v2_vectet = vtetat2 - vtetat1
            v2 = m.sqrt(v2_vecr**2+v2_vectet**2)  #Magnitude of Delta V2

            #S/C velocities at perigee on transfer orbit 2
            #nup = m.acos(1/e*(p/rp-1))
            nupt2 = 0
            nupt2_deg = nupt2*180/m.pi
            vrpt2 = self.mu*et2*m.sin(nupt2)/ht2
            vtetpt2 = ht2/rpt2

            #S/C velocities on final orbit
            self.vr2 = 0  #Spacecraft velocity is entirely in theta hat direction
            vtet2 = self.h2/self.r2 
            v2 = m.sqrt(self.mu2/self.r2) #Velocity magnitude anywhere on the final circular orbit

            #Delta V3
            v3_vecr = self.vr2 - vrpt2  #Impulsive burn from transfer orbit to final orbit
            v3_vectet = vtet2 - vtetpt2
            v3 = m.sqrt(v3_vecr**2+v3_vectet**2)  #Magnitude of Delta V2

            #Calculating Time of Flight in MONTHS AND DAYS for Transfer 1
            TOFmd1 = m.pi*m.sqrt(at1**3/self.mu)/3600/24/30 #Approximately calculates time of flight in months
            TOFm1 = int(TOFmd1//1) #Calculates TOF to nearest month
            TOFr11 = TOFmd1%1  #Calculates remainder of above TOF calculation
            TOFd1 = int(TOFr11*30) #Finds number of days left in TOF calculation
            #Calculating Time of Flight in MONTHS AND DAYS for Transfer 2
            TOFmd2 = m.pi*m.sqrt(at2**3/self.mu)/3600/24/30 #Approximately calculates time of flight in months
            TOFm2 = int(TOFmd2//1) #Calculates TOF to nearest month
            TOFr12 = TOFmd2%1  #Calculates remainder of above TOF calculation
            TOFd2 = int(TOFr12*30) #Finds number of days left in TOF calculation
                
            #Calculating Time of Flight in HOURS AND MINUTES for Transfer 1
            TOFhm1 = m.pi*m.sqrt(at1**3/self.mu)/3600  #Exactly calculates time of flight in hours
            TOFh1 = int(TOFhm1//1)  #Calculates TOF to nearest hour
            TOFr21 = TOFhm1%1 #Calculates remainder of above TOF calculation
            TOFmn1 = int(TOFr21*60)  #Finds number of minutes left in TOF calculation
            #Calculating Time of Flight in HOURS AND MINUTES for Transfer 2
            TOFhm2 = m.pi*m.sqrt(at2**3/self.mu)/3600  #Exactly calculates time of flight in hours
            TOFh2 = int(TOFhm2//1)  #Calculates TOF to nearest hour
            TOFr22 = TOFhm2%1 #Calculates remainder of above TOF calculation
            TOFmn2 = int(TOFr22*60)  #Finds number of minutes left in TOF calculation

            #Calculating Total Time of Flight in MONTHS AND DAYS
            TOFmT = TOFm1 + TOFm2
            TOFdT = TOFd1 + TOFd2
            #Calculating Total Time of Flight in HOURS AND MINUTES
            TOFhT = TOFh1 + TOFh2
            TOFmnT = TOFmn1 + TOFmn2

            if TOFmT == 0 and TOFdT < 4:
                TOF1 = TOFhT
                TOF2 = TOFmnT
                unit1 = 'hours'
                unit2 = 'minutes'
            else:
                TOF1 = TOFmT
                TOF2 = TOFdT
                unit1 = 'months'
                unit2 = 'days'

            #Calculating mass of propellant needed for mission (all in kg)
            mf3 = self.m_sc  #The spacecraft's final mass after burn 2 (equal to dry mass)
            mp3 = mf3*(m.e**(v3/(self.isp*self.g))-1)
            mf2 = mp3 + mf3
            mp2 = mf2*(m.e**(v2/(self.isp*self.g))-1) #Propellant mass needed to execute burn 2
            mf1 = mp2 + mf2  #Propellant mass needed after burn 1
            mp1 = mf1*(m.e**(v1/(self.isp*self.g))-1) #Propellant mass need to execute burn 1
            mpT = mp1 + mp2 + mp3  #Total propellant mass need for whole mission
            mo = mpT + self.m_sc  #Total mass of S/C at beginning of mission (throw mass)
            MF = mf3/mo*100  #Mass Fraction as a percentage

            #Printing results for Bi-Elliptic Transfer
            if  MF > 10:
                ellip_opt_1 = "The Bi-Elliptic transfer" 
                ellip_opt_2 = "IS optimal for this transfer."
            else:
                ellip_opt_1 = "The Bi-Elliptic transfer is" 
                ellip_opt_2 = "NOT optimal for this transfer."

            self.print_ellip_opt_1 = ttk.Label(self.mainframe, text=ellip_opt_1, foreground='white', background='black', font=("Brass Mono", 15))
            self.print_ellip_opt_1.grid(row=10, column=0, pady=15)
            self.print_ellip_opt_2 = ttk.Label(self.mainframe, text=ellip_opt_2, foreground='white', background='black', font=("Brass Mono", 15))
            self.print_ellip_opt_2.grid(row=10, column=1, pady=15)

            #Using a more reader-friendly value for semi-major axis of both half-orbits and mass of propellant
            at1_rep = at1  #Reported semi-major axis for transfer orbit 1
            at2_rep = at2  #Reported semi-major axis for transfer orbit 2
            mpT_rep = mpT  #Reported mass of propellant required
            exp_at1 = 0  #Reported exponent for semi-major axis measurement 1
            exp_at2 = 0  #Reported exponent for semi-major axis measurement 2
            exp_mp = 0  #Reported exponent for mass propellant measurement
            for n in range(0,25):
                if not at1_rep < 10:
                    at1_rep = at1/(10**n)
                    exp_at1 = n
                if not at2_rep < 10:
                    at2_rep = at2/(10**n)
                    exp_at2 = n
                if not mpT_rep < 10:
                    mpT_rep = mpT/(10**n)
                    exp_mp = n
                    
            ######################
            # Labels for Results #
            ######################
            
            at1_text1 = "Semi-Major Axis (a)\nfor Transfer Orbit 1:"
            at1_text2 = f"{round(at1_rep,3)}*10^{exp_at1} km"
            self.print_at1_1 = ttk.Label(self.mainframe, text=at1_text1, foreground='white', background='black', font=("Brass Mono", 14))
            self.print_at1_1.grid(row=11, column=0, pady=1)
            self.print_at1_2 = ttk.Label(self.mainframe, text=at1_text2, foreground='white', background='black', font=("Brass Mono", 15))
            self.print_at1_2.grid(row=12, column=0, pady=1)

            et1_text1 = "Eccentricity (e)\nfor Transfer Orbit 1:"
            et1_text2 = f"{round(et1,3)}"
            self.print_et1_1 = ttk.Label(self.mainframe, text=et1_text1, foreground='white', background='black', font=("Brass Mono", 14))
            self.print_et1_1.grid(row=11, column=1, pady=1)
            self.print_et1_2 = ttk.Label(self.mainframe, text=et1_text2, foreground='white', background='black', font=("Brass Mono", 15))
            self.print_et1_2.grid(row=12, column=1, pady=1)

            at2_text1 = "Semi-Major Axis (a)\nfor Transfer Orbit 2:"
            at2_text2 = f"{round(at2_rep,3)}*10^{exp_at2} km"
            self.print_at2_1 = ttk.Label(self.mainframe, text=at2_text1, foreground='white', background='black', font=("Brass Mono", 14))
            self.print_at2_1.grid(row=11, column=2, pady=1)
            self.print_at2_2 = ttk.Label(self.mainframe, text=at2_text2, foreground='white', background='black', font=("Brass Mono", 15))
            self.print_at2_2.grid(row=12, column=2, pady=1)

            et2_text1 = "Eccentricity (e)\nfor Transfer Orbit 2:"
            et2_text2 = f"{round(et2,3)}"
            if et2 < 0:
                et2_text2 = f"{round(-et2,3)}"
            self.print_et2_1 = ttk.Label(self.mainframe, text=et2_text1, foreground='white', background='black', font=("Brass Mono", 14))
            self.print_et2_1.grid(row=11, column=3, pady=1)
            self.print_et2_2 = ttk.Label(self.mainframe, text=et2_text2, foreground='white', background='black', font=("Brass Mono", 15))
            self.print_et2_2.grid(row=12, column=3, pady=1)
            
            deltv1_text1 = "ΔV1:"
            deltv1_text2 = f"{round(v1,2)} km/s"
            self.print_deltv1_1 = ttk.Label(self.mainframe, text=deltv1_text1, foreground='white', background='black', font=("Brass Mono", 15))
            self.print_deltv1_1.grid(row=13, column=1, pady=15)
            self.print_deltv1_2 = ttk.Label(self.mainframe, text=deltv1_text2, foreground='white', background='black', font=("Brass Mono", 15))
            self.print_deltv1_2.grid(row=14, column=1, pady=1)

            deltv2_text1 = "ΔV2:"
            deltv2_text2 = f"{round(v2,2)} km/s"
            self.print_deltv2_1 = ttk.Label(self.mainframe, text=deltv2_text1, foreground='white', background='black', font=("Brass Mono", 15))
            self.print_deltv2_1.grid(row=13, column=2, pady=15)
            self.print_deltv2_2 = ttk.Label(self.mainframe, text=deltv2_text2, foreground='white', background='black', font=("Brass Mono", 15))
            self.print_deltv2_2.grid(row=14, column=2, pady=1)

            deltv3_text1 = "ΔV3:"
            deltv3_text2 = f"{round(v3,2)} km/s"
            self.print_deltv3_1 = ttk.Label(self.mainframe, text=deltv3_text1, foreground='white', background='black', font=("Brass Mono", 15))
            self.print_deltv3_1.grid(row=13, column=3, pady=15)
            self.print_deltv3_2 = ttk.Label(self.mainframe, text=deltv3_text2, foreground='white', background='black', font=("Brass Mono", 15))
            self.print_deltv3_2.grid(row=14, column=3, pady=1)
            
            TOF_text1 = "Time of Flight (TOF):"
            TOF_text2 = f"{round(TOF1,2)} {unit1}\n{round(TOF2,2)} {unit2}"
            self.print_TOF_1 = ttk.Label(self.mainframe, text=TOF_text1, foreground='white', background='black', font=("Brass Mono", 14))
            self.print_TOF_1.grid(row=13, column=0, pady=1)
            self.print_TOF_2 = ttk.Label(self.mainframe, text=TOF_text2, foreground='white', background='black', font=("Brass Mono", 15))
            self.print_TOF_2.grid(row=14, column=0, pady=1)

            self.sc_label = Label(self.mainframe, text=f'Spacecraft Used:', foreground='White', background='Black', font=("Brass Mono", 17))
            self.sc_label.grid(row=15,column=0,sticky='NWES')
            
            self.sc_img = Image.open(f"2d_transfer_Pics/{self.set_spacecraft.get()}_Pic.png")
            self.sc_img_resized = self.sc_img.resize((100,100), Image.ANTIALIAS)
            self.new_sc_img = ImageTk.PhotoImage(self.sc_img_resized)
            self.sc_img_label = Label(self.mainframe, image=self.new_sc_img, background='Black')
            self.sc_img_label.grid(row=15, column=1, sticky='NWES')

            self.sc_name = Label(self.mainframe, text=f'{self.set_spacecraft.get()}', foreground='White', background='Black', font=("Brass Mono", 15))
            self.sc_name.grid(row=16,column=1,sticky='NWES')
            
            mpT_text1 = "Propellant Required:"
            mpT_text2 = f"{round(mpT_rep,3)}*10^{exp_mp} kg"
            self.print_mpT_text1 = ttk.Label(self.mainframe, text=mpT_text1, foreground='white', background='black', font=("Brass Mono", 14))
            self.print_mpT_text1.grid(row=15, column=2, pady=1)
            self.print_mpT_text2 = ttk.Label(self.mainframe, text=mpT_text2, foreground='white', background='black', font=("Brass Mono", 15))
            self.print_mpT_text2.grid(row=16, column=2, pady=1)

            MF_text1 = "Mass Fraction (MF):"
            MF_text2 = f"{round(MF,2)}%"
            self.print_MF_text1 = ttk.Label(self.mainframe, text=MF_text1, foreground='white', background='black', font=("Brass Mono", 14))
            self.print_MF_text1.grid(row=15, column=3, pady=15)
            self.print_MF_text2 = ttk.Label(self.mainframe, text=MF_text2, foreground='white', background='black', font=("Brass Mono", 15))
            self.print_MF_text2.grid(row=16, column=3, pady=1)
            

            #Plotting the spacecraft's bi-elliptic transfer orbit
            t = np.linspace(-2*m.pi,2*m.pi,10**6)
            z = np.linspace(0,m.pi,10**6)
            zan = np.linspace(0,np.pi,100)  #Sets values for animation functionto iterate thru
            zz = np.linspace(m.pi,2*m.pi,10**6)
            zzan = np.linspace(-np.pi,0,100)

            #Transfer orbit 1
            if self.r1 < self.r2:
                xt1 = (at1*et1+rpt1)*np.cos(z)-1/2*(rat1-rpt1)
                xt1an = (at1*et1+rpt1)*np.cos(zan)-1/2*(rat1-rpt1)  #Sets values for animation function to iterate thru
            elif self.r1 > self.r2:
                xt1 = (at1*et1+rpt1)*np.cos(z)-1/2*(rat1-rpt1)
                xt1an = (at1*et1+rpt1)*np.cos(zan)-1/2*(rat1-rpt1)  #Sets values for animation function to iterate thru
            yt1 = bt1*np.sin(z) 
            yt1an = bt1*np.sin(zan)  #Sets values for animation function to iterate thru
            #Transfer orbit 2
            if self.r1 < self.r2:
                xt2 = (at2*et2+rpt2)*np.cos(zz)-1/2*(rat2-rpt2)
                xt2an = (at2*et2+rpt2)*np.cos(zzan)-1/2*(rat2-rpt2)  #Sets values for animation function to iterate thru
            elif self.r1 > self.r2:
                xt2 = (at2*et2+rpt2)*np.cos(zz)-1/2*(rat2-rpt2)
                xt2an = (at2*et2+rpt2)*np.cos(zzan)-1/2*(rat2-rpt2)  #Sets values for animation function to iterate thru
            yt2 = bt2*np.sin(zz) 
            yt2an = bt2*np.sin(zzan)  #Sets values for animation function to iterate thru
            mark = 8

            plt.style.use('dark_background')
            fig, ax = plt.subplots(figsize=(9,9))
            #ax = plt.axes()
            plt.xlabel("Distance (km)")
            plt.ylabel("Distance (km)")
            
            if self.trsf_type == '1':
                plt.title(f"{self.sc_type} Spacecraft's Bi-Elliptic Transfer \nfrom {self.planet} Orbit at Altitude of {round(self.alt1,2)} km to Altitude of {round(self.alt2,2)} km")
                #Central Body
                ax.plot(0,0,marker='o',
                    markersize = 25, color = self.clr)

            elif self.trsf_type == '2':
                plt.title(f"{self.sc_type} Spacecraft's Bi-Elliptic Transfer \nfrom {self.planet1} to {self.planet2}")
                #Sun
                ax.plot(0,0, marker='o',
                    markersize = 25, color='darkorange')
            elif self.trsf_type == 'pm' or self.trsf_type == 'mp':
                plt.title(f"{self.sc_type} Spacecraft's Bi-Elliptic Transfer \nfrom {self.body1} to {self.body2}")
                #Planet
                if self.trsf_type == 'pm':
                    color = self.clr1
                elif self.trsf_type == 'mp':
                    color = self.clr2
                ax.plot(0,0, marker='o',
                    markersize = 25, color=color)
            
            #Orbit around Planet 1
            x1 = self.r1*np.cos(t)
            y1 = self.r1*np.sin(t)
            ax.plot(x1,y1, color='lightsteelblue')
            #Orbit around Planet 2
            x2 = self.r2*np.cos(t)
            y2 = self.r2*np.sin(t)
            ax.plot(x2,y2, color='lightsteelblue')
            #Transfer Orbit 1 & 2
            ax.plot(xt1,yt1,label = "Transfer Orbit 1",color = 'red')
            ax.plot(xt2,yt2,label = "Transfer Orbit 2",color = 'lime')
            
            #Spacecraft Locations 1 and 2
            if self.trsf_type == '1':
                ax.plot(-self.r1*np.sin(nu1),0, marker=(4,0,0),
                    markersize = 12, color='yellow',
                        label = 'Initial S/C Position', alpha=0.4)
                ax.plot(-self.r2*np.sin(nu2),0, marker=(4,0,0),
                    markersize = 12, color='lime',
                        label = 'Final S/C Position', alpha=0.4)
            elif self.trsf_type == '2':  #This plots the location of each body on their respective orbtis
                ax.plot(-self.r1*np.sin(nu1),0, marker='o',
                    markersize = 6.378*3, color=self.clr1, label = f'{self.planet1}')
                ax.plot(-self.r2*np.sin(nu2),0, marker='o',
                    markersize = 6.378*3, color=self.clr2, label = f'{self.planet2}')
            elif self.trsf_type == 'pm':
                ax.plot(-self.r1*np.sin(nu1),0, marker=(4,0,0),
                    markersize = 12, color='yellow',
                        label = 'Initial S/C Position', alpha=0.4)
                ax.plot(-self.r2*np.sin(nu2),0, marker='o',
                    markersize = 6.378*3, color=self.clr2, label = self.body2)
            elif self.trsf_type == 'mp':
                ax.plot(-self.r1*np.sin(nu1),0, marker='o',
                    markersize = 6.378*3, color=self.clr1,
                        label = self.body1)
                ax.plot(-self.r2*np.sin(nu2),0, marker=(4,0,0),
                    markersize = 12, color='lime',
                        label = 'Final S/C Position', alpha=0.4)

            #Plotting spacecraft's movement along transfer orbit
            trsf1, = ax.plot([],[],marker=(4,0,0),
                            color=self.scclr,markersize=12,label=f'{self.sc_type} Spacecraft')
            trsf2, = ax.plot([],[],marker=(4,0,0),
                            color=self.scclr,markersize=12)
            def animate1(i):
                trsf1.set_data(xt1an[i], yt1an[i])
                return trsf1
            def animate2(i):
                trsf2.set_data(xt2an[i], yt2an[i])
                return trsf2

            frames = 100
            interval = 200
            plt.grid(color = 'lightgray', linestyle='--',linewidth = 0.3)
            plt.legend(bbox_to_anchor=(0.95, 0.85), fontsize='16.5', loc='lower center')
            plt.pause(3)
            anim1 = FuncAnimation(fig, animate1, frames=frames, interval=interval, repeat=True, repeat_delay = frames*interval)
            plt.pause((frames+2)*interval/10**3)
            anim2 = FuncAnimation(fig, animate2, frames=frames, interval=interval, repeat=False, repeat_delay = (frames+2)*interval)
            plt.pause((frames+2)*interval/10**3)
            plt.show()
            
            
        else:
            print("¯\_(ツ)_/¯")

warnings.filterwarnings("ignore")
    
if __name__ == '__main__':
    TransferSim()
