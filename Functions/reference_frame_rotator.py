# ------------------------ Reference Frame Rotations ---------------------------- #

import numpy as np
import matplotlib.pyplot as plt
import spiceypy as spice

spice.furnsh('spice_meta_kernal.txt')
spice.furnsh('earth_000101_240603_240310.bpc')

def choose_DCM(DCM, theta):

    DCMs = {
            # Roll Rotation
            'x' : np.matrix(
                    [[1, 0, 0], 
                    [0, np.cos(theta), -np.sin(theta)],
                    [0, np.sin(theta), np.cos(theta)]]
                ),

            # Pitch Rotation
            'y' : np.matrix(
                    [[np.cos(theta), 0, np.sin(theta)], 
                    [0, 1, 0],
                    [-np.sin(theta), 0, np.cos(theta)]]
                ),

            # Yaw Rotation
            'z' : np.matrix(
                    [[np.cos(theta), -np.sin(theta), 0], 
                    [np.sin(theta), np.cos(theta), 0],
                    [0, 0, 1]]
                ),
        }

    return DCMs[f'{DCM}']

def eci2ecef(rs, tspan, frame= 'J2000'):

    ns = rs.shape[0]
    Cs = np.zeros((ns, 3, 3))       # Holds all rotation matrices
    rs_ecef = np.zeros(rs.shape)    # Holds ECEF position vectors

    for n in range(ns):
        Cs[n, :, :] = spice.pxform(frame, 'ITRF93', tspan[n])   # ITRF93 is how ECEF frame is defined in Spice
        rs_ecef[n, :] = np.dot(Cs[n, :, :], rs[n, :])           # This is where you multiply your rotation matrix and ECI vector for each step n

    return rs_ecef, Cs

def eci2latlong(rs, tspan, frame= 'J2000'):
    
    ns = rs.shape[0]
    lat_longs = np.zeros((ns, 3))
    rs_ecef, Cs = eci2ecef(rs, tspan)

    for n in range(ns):
        r_mag, long, lat = spice.reclat(rs_ecef[n, :])
        lat_longs[n,:] = [np.degrees(lat), np.degrees(long), r_mag]

    return lat_longs, rs_ecef, Cs


class frame_rotation():

    def __init__(self, original_vector,**kwargs):
        
        self.original_vector = original_vector.reshape((3,1))
        
        theta = 0
        if 'theta' in kwargs:
            theta = kwargs['theta']
        self.theta = np.radians(theta)

        self.principle_axis = 'z'
        if 'principle' in kwargs:
            self.principle_axis = kwargs['principle']

        self.thetas = []
        if 'thetas' in kwargs:
            theta1 = kwargs['thetas'][0]   
            theta2 = kwargs['thetas'][1]
            theta3 = kwargs['thetas'][2]
            self.thetas = [np.radians(theta1), np.radians(theta2), np.radians(theta3)]

        self.spice = True
        if self.spice == False:
            pass
        self.spiceDCM = []

        self.euler = 0
        if 'euler' in kwargs:
            self.euler = kwargs['euler']


    def find_rotation_matrix(self, **kwargs):
        
        RAAN = np.radians(30)   # RAAN for a Molniya Orbit
        inc = np.radians(63.4)  # Inclination for Molniya Orbit
        w = np.radians(270)     # Argument of Periapsis for Molniya Orbit

        euler = self.euler

        if not self.thetas == []:
            if 'thetas' in kwargs:
                RAAN = np.radians(kwargs['thetas'][0])
                inc = np.radians(kwargs['thetas'][1])
                w = np.radians(kwargs['thetas'][2])
            else:
                RAAN, inc, w = self.thetas[0], self.thetas[1], self.thetas[2]

        if 'euler' in kwargs:
            euler = kwargs['euler']

        first = ''
        second = ''
        third = ''
        if euler == 321:
            first = 1
            second = 2
            third = 3
        if euler == 313:
            first = 3
            second = 1
            third = 3

        RM = spice.eul2m(w, inc, RAAN, first, second, third)   # RM = Rotation Matrix; Order for this function is third rotation, then second, then first (euler 321 would be 1, 2, 3)
        self.spiceDCM = RM.T

        return self.spiceDCM


    def rotate_vector(self, **kwargs):
        
        self.DCM = choose_DCM(self.principle_axis, self.theta)
        if 'spice' in kwargs:
            self.DCM = self.spiceDCM

        self.new_vector = self.DCM @ self.original_vector

        self.new_x = 0
        self.new_y = 0
        self.new_z = 0
        if not self.euler == 0:
            self.new_vector = self.Euler()

            self.new_x = self.Euler(vector= np.array([1,0,0]))
            #self.new_x[0] = -1 * self.new_x[0]      # For plotting; keeps x axis in correct direction
            self.new_y = self.Euler(vector= np.array([0,1,0]))
            #self.new_y[1] = -1 * self.new_y[1]
            self.new_z = self.Euler(vector= np.array([0,0,1]))
            #self.new_z[2] = -1 * self.new_z[2]

        return self.new_vector

    
    def plot_new_frame(self):

        spiceDCM = self.find_rotation_matrix(euler= self.euler, thetas= self.thetas)

        plt.style.use('dark_background')
        ax = plt.figure().add_subplot(projection= '3d')

        l = 100
        x, y, z = [[l,0,0], [0,l,0], [0,0,l]]
        zeros = [[0,0,0], [0,0,0], [0,0,0]]
        xprime, yprime, zprime = spiceDCM[0, :] * l, spiceDCM[1, :] * l, spiceDCM[2, :] * l

        # PLOTTING AXES #
        ax.quiver(zeros, zeros, zeros, x, y, z , color= 'crimson', label= 'Original Coordinate Frame')    # Original x-y-z axes
        ax.text( 1.1*l, 0, 0, 'X', color = 'crimson' )
        ax.text( 0, 1.1*l, 0, 'Y', color = 'crimson'  )
        ax.text( 0, 0, 1.1*l, 'Z', color = 'crimson' )
        ax.quiver(zeros, zeros, zeros, xprime, yprime, zprime, color= 'lime', alpha= 0.25, label= 'New Coordinate Frame')
        ax.text( spiceDCM[ 0, 0 ] * 1.1*l + 0.2*l, spiceDCM[ 1, 0 ] * 1.1*l, spiceDCM[ 2, 0 ] * 1.1*l, " X' ",
				color = 'lime' )
        ax.text( spiceDCM[ 0, 1 ] * 1.1*l, spiceDCM[ 1, 1 ] * 1.1*l  + 0.2*l, spiceDCM[ 2, 1 ] * 1.1*l, " Y' ",
            color = 'lime' )
        ax.text( spiceDCM[ 0, 2 ] * 1.1*l, spiceDCM[ 1, 2 ] * 1.1*l, spiceDCM[ 2, 2 ] * 1.1*l  + 0.2*l, " Z' ",
            color = 'lime' )

        plt.xlim(-l*1.35, l*1.35)
        plt.ylim(-l*1.35, l*1.35)
        ax.set_zlim(-l*1.35, l*1.35)
        
        plt.legend()
        plt.show()
