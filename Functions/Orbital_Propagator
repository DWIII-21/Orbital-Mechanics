# --------------------- Orbit Propagator -------------------- #

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import planetary_data as pd
import class2cartcomponents as cl2crt

'''
The purpose of this script is to create an RK4 integrator to gain a better understanding of how they work.
It is *almost* always optimal to use "solve_ivp" for an practical applications.
'''

def rk4singlestep(function, dt, t0, y0):

    '''
    This function is a single 4th order Runge-Kutta integrator where:
    function = the ODE in question
    self.dt = time step
    t0 = initial t condition
    y0 = initial y condition as a state vector (contains xdot, ydot, and zdot); also written as y_k
    '''

    k1 = function(t0, y0)
    k2 = function(t0 + dt/2, y0 + (dt/2)*k1)
    k3 = function(t0 + dt/2, y0 + (dt/2)*k2)
    k4 = function(t0 + dt, y0 + dt*k3)

    y1 = y0 + (dt/6)*(k1 + 2*k2 + 2*k3 + k4)  # this is y_k+1

    return y1


def get_bodies(body):
    
    # Available bodies
    if body == 'sun':
        body = pd.sun
    if body == 'mercury':
        body = pd.mercury
    if body == 'venus':
        body = pd.venus
    if body == 'earth':
        body = pd.earth
    if body == 'moon':
        body = pd.moon
    if body == 'mars':
        body = pd.mars
    if body == 'jupiter':
        body = pd.jupiter
    if body == 'saturn':    
        body = pd.saturn
    if body == 'titan':    
        body = pd.titan
    if body == 'uranus':    
        body = pd.uranus
    if body == 'neptune':    
        body = pd.neptune
    if body == 'pluto':    
        body = pd.pluto

    return body

def CR3BP_system(system_name):

        SYSTEMS = {
                'earth-moon' : {
                    'mu' : pd.moon['mu'] / ( pd.earth['mu'] + pd.moon['mu']),         # m2/(m1+m2)  (NOT gravitational parameter)
                    'L1' : 0.8362925909457339,
                    'L2' : 1.1561681659055247,
                    'L3' : -1.0051155116068917,
                    'body1' : pd.earth,
                    'body2' : pd.moon,
                    'r1' : 4671,   # Distance from Earth to the Earth-Moon Barycenter in km
                    'r2' : 382138.98904109595,  # Distance from Moon to Earth-Moon Barycenter in km
                },

                'sun-jupiter': {
                    'mu' : 0.000953875
                }
            }

        return SYSTEMS[f'{system_name}']

#############################
# 2-BODY RK4 ODE PROPAGATOR #
#############################
  
class rk4_2Body():

    '''
    Uses an RK4 ODE solver to propagate and plot orbits in the Inertial Reference Frame.
    '''

    def __init__(self, **kwargs):

        '''
        Establishes parameters necessary for calculations:
        cb = cb that is being orbited
        alt = altitude above cb at which object is orbiting
        '''

        # Setting Default Parameters

        alt = 1000          # default altitude in meters
        self.mu = 1
        self.eq = 100
        self.gravity = 5
        self.color = 'silver'
        self.map_color = 'binary'
        self.name = 'Central Body'
        self.big_body = False
        self.show = True

        self.dt = 100
        self.tspan = 100 * 60  # Upper bound of integration (from 0 to T)

        # Establishing "cb" keyword argument

        if 'cb' in kwargs:      # Defines the central body of the system_name
            cb = get_bodies(kwargs['cb'])          

            self.mu = cb['mu']
            self.eq = cb['eq']
            self.gravity = cb['grav']
            self.color = cb['clr']
            self.map_color = cb['clrmap']
            self.name = cb['name']

        # Establishing "Altitude" and "Time" keyword arguments
            
        if 'alt' in kwargs:
            alt = kwargs['alt']

        if 'dt' in kwargs:
            self.dt = kwargs['dt']

        if 'tspan' in kwargs:
            self.tspan = kwargs['tspan']


        self.alt = alt
        self.h = self.eq + self.alt

        r_mag = self.h
        v_mag = np.sqrt(self.mu / self.h)
        # Set initial conditions
        y0 = [r_mag, 0, 0, 0, v_mag, 0]

        if 'classical' in kwargs:       # This argument must be entered as a tuple
            a = kwargs['classical'][0]
            e = kwargs['classical'][1]
            i = kwargs['classical'][2]
            omega = kwargs['classical'][3]
            w = kwargs['classical'][4]
            M = kwargs['classical'][5]

            x, y, z, xdot, ydot, zdot = cl2crt.class2cartcomp(self.mu, a, e, i, omega, w, M)

            y0 = [x, y, z, xdot, ydot, zdot]

        elif 'y0' in kwargs:    # Takes input of state vector in inertial reference frame
            y0 = kwargs['y0']

        self.y0 = y0

        self.max_step = np.inf
        if 'maxstep' in kwargs:
            self.max_step = kwargs['maxstep']

        # Establishing keyword arguments associated with plotting 

        if 'bb' in kwargs:
            self.big_body = kwargs['bb']

        if 'show' in kwargs:
            self.show = kwargs['show']

        self.object_label = 'Spacecraft'
        if 'objlabel' in kwargs:
            self.object_label = kwargs['objlabel']

        self.start_label = f'Initial Position of {self.object_label}'
        if 'startpoint' in kwargs:
            self.start_label = kwargs['startpoint']

        self.title = f'Trajectory of {self.object_label}'
        if 'title' in kwargs:
            self.title = kwargs['title']

        self.axes = False
        if 'axes' in kwargs:
            self.axes = kwargs['axes']
    

    def state_diff_eq(self, t, y):
        
        rx, ry, rz, vx, vy, vz = y

        r_vector = np.array([rx, ry, rz])
        r_magnitude = np.linalg.norm(r_vector)

        ax, ay, az = -(r_vector * self.mu / r_magnitude**3)     # From EOM: ȑ_vector + (μ/r^3)*r_vector = 0

        return [vx, vy, vz, ax, ay, az]
    
    def propagate_orbit(self):
        
        print(f'Propagating orbit...')

        # Compute trajectory
        number_of_steps = int(self.tspan/self.dt)
        t = np.linspace(0, self.tspan, number_of_steps)

        # Plotting results from one of Python's built in RK4 solver "solve_ivp"
        orbit_solution = solve_ivp(
            fun = self.state_diff_eq, 
            t_span = (0, self.tspan), 
            y0 = self.y0, 
            t_eval = t,
            max_step = self.max_step,
            rtol = 3e-14,
            atol = 3e-14,
        )
        t = orbit_solution.t
        y = orbit_solution.y.T

        if self.show == False:
            return t, y
        else:
            # Plotting results from custom RK4 solver
            plt.style.use('dark_background')
            ax = plt.figure().add_subplot(projection= '3d')

            # Plotting Central Body
            if self.big_body:
                # Plotting Large Spherical Central Body
                _u,_v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
                _x = self.eq * np.cos(_u)*np.sin(_v)
                _y = self.eq * np.sin(_u)*np.sin(_v)
                _z = self.eq * np.cos(_v)

                ax.plot_surface(_x,_y,_z, cmap= self.map_color)
            else:
                # Plotting Central Body as a Point Mass (Dot)
                ax.plot(0, 0, 0, color= self.color, marker= '.', markersize= 21, label= self.name)

            # Plotting Orbit about Central Body
            ax.plot(y[:, 0], y[:, 1], y[:, 2], color= 'orchid', label= f'Orbit of {self.object_label} about {self.name}')
            ax.plot(self.y0[0], self.y0[1], self.y0[2], color= 'magenta', marker= '.', markersize= 10, label= self.start_label)

            if self.axes:
                # Plot x,y,z axis vectors
                l = abs(self.y0[0] * 0.5)
                x,y,z = [[0,0,0], [0,0,0], [0,0,0]]
                u,v,w = [[l,0,0], [0,l,0], [0,0,l]]

                ax.quiver(x,y,z,u,v,w, color = 'lime')

            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_zlabel('Z')
            plt.legend()
            plt.show()

###############
# CR3BP CLASS #
###############

class CR3BP():

    '''
    Propagates and plots orbits under restrictions and assumptions from the Circular Restricted Three Body Problem.
    Plots are done with respect to the Rotating Reference Frame.
    '''
    
    def __init__(self, **kwargs):

        # Setting Default Parameters

        self.eq = 100
        self.gravity = 5
        self.color = 'silver'
        self.map_color = 'binary'
        self.name = 'Central Body'
        self.big_body = False
        self.plot = False

        self.b1 = pd.sun
        self.b2 = pd.earth

        self.dt = 100
        self.tspan = 100 * 60  # Upper bound of integration (from 0 to T)

        if 'system' in kwargs:
            system_name = kwargs['system']
        
            if system_name == 'em' or system_name == 'earth' or system_name == 'moon' or system_name == 'earth-moon':
                system_name = 'earth-moon'
                self.b1 = pd.earth
                self.b2 = pd.moon
                r1 = CR3BP_system(system_name)['r1']
                r2 = CR3BP_system(system_name)['r2']
            if system_name == 'sj' or system_name == 'jupiter':
                system_name = 'sun-jupiter'
                self.b1 = pd.sun
                self.b2 = pd.jupiter
            
            self.sys = CR3BP_system(system_name)
            self.mu = self.sys['mu']

        if 'body1' in kwargs:
            self.b1 = get_bodies(kwargs['body1'])
        if 'body2' in kwargs:
            self.b2 = get_bodies(kwargs['body2'])
        
        self.mu = self.b2['mass'] / (self.b1['mass'] + self.b2['mass'])

        if 'dt' in kwargs:
            self.dt = kwargs['dt']

        if 'tspan' in kwargs:
            self.tspan = kwargs['tspan']

        if 'mu' in kwargs:
            self.mu = kwargs['mu']
        
        y0 = []
        if 'classical' in kwargs:       # This argument must be entered as a tuple
            mu_gravparam = kwargs['classical'][0]   # This value must be the Gravitational Parameter of the Largest/Central Body
            a = kwargs['classical'][1]
            e = kwargs['classical'][2]
            i = kwargs['classical'][3]
            omega = kwargs['classical'][4]
            w = kwargs['classical'][5]
            M = kwargs['classical'][6]

            x, y, z, xdot, ydot, zdot = cl2crt.class2cartcomp(mu_gravparam, a, e, i, omega, w, M)

            y0 = [x, y, z, xdot, ydot, zdot]

        elif 'y0' in kwargs:    # Takes input of initial state vector in rotating reference frame (radii must be relative to mu of system, velocity values unchanged from inertial frame)
            y0 = kwargs['y0']

        elif 'y0km' in  kwargs:     # Takes input of initial state vector in inertial reference frame (radii in km)
            y0 = kwargs['y0km']
            x, y, z, xdot, ydot, zdot = y0
            y0 = [x / (r1+r2), y / (r1+r2), z / (r1+r2), xdot, ydot, zdot]

        self.y0 = y0

        self.max_step = np.inf
        if 'maxstep' in kwargs:
            self.max_step = kwargs['maxstep']


        # Establishing keyword arguments for specific position pertabations (used in determining orbit stability)
        self.dx = 0
        self.dy = 0
        self.dz = 0
        if 'dx' in kwargs:
            self.dx = kwargs['dx']
        if 'dy' in kwargs:
            self.dy = kwargs['dy']
        if 'dz' in kwargs:
            self.dz = kwargs['dz']

        # Establishing keyword arguments associated with plotting

        if 'plot' in kwargs:
            self.plot = kwargs['plot']

        self.object_label = 'Spacecraft'
        if 'objlabel' in kwargs:
            self.object_label = kwargs['objlabel']

        self.start_label = f'Initial Position of {self.object_label}'
        if 'startpoint' in kwargs:
            self.start_label = kwargs['startpoint']

        self.title = f'Trajectory of {self.object_label}'
        if 'title' in kwargs:
            self.title = kwargs['title']
    
    def state_diff_eq(self, t, y):

        rx, ry, rz, vx, vy, vz = y
        rx += self.dx
        ry += self.dy
        rz += self.dz

        r13_vec = [rx + self.mu, ry, rz]
        r13_mag = np.linalg.norm(r13_vec)       # Magnitude of R_13 vector
        r23_vec = [rx - 1 + self.mu, ry, rz]
        r23_mag = np.linalg.norm(r23_vec)

        # Gradient of the Pseudo-Potential Function
        omega_x = rx  -  (1 - self.mu) * (rx + self.mu) / r13_mag**3  -  self.mu * (rx - 1 + self.mu) / r23_mag**3
        omega_y = ry  -  (1 -self.mu) * ry / r13_mag**3  -  self.mu * ry / r23_mag**3
        omega_z = -(1 - self.mu) * rz / r13_mag**3  -  self.mu * rz / r23_mag**3

        ax = 2 * vy + omega_x
        ay = -2 * vx + omega_y
        az = omega_z

        return [vx, vy, vz, ax, ay, az]
    
    def propagate_orbit(self):

        print(f'Propagating orbit...')

        # Compute trajectory
        number_of_steps = int(self.tspan/self.dt)
        t = np.linspace(0, self.tspan, number_of_steps)

        # Plotting results from one of Python's built in RK4 solver "solve_ivp" using the "LSODA" method
        orbit_solution = solve_ivp(
            fun = self.state_diff_eq, 
            t_span = (0, self.tspan), 
            y0 = self.y0, 
            t_eval = t,
            method = 'LSODA',
            max_step = self.max_step,
            rtol = 3e-14,
            atol = 3e-14, )
        
        ts = orbit_solution.t
        ys = orbit_solution.y.T
        n_steps = orbit_solution.y.shape[0]

        if self.plot == '2d':
            # Plotting results for CR3BP
            plt.style.use('dark_background')
            plt.figure()
            plt.plot(-self.mu, 0, marker= '.', markersize= 18, color= self.b1['clr'], label= self.b1['name'])     # Plotting Large Body (Body 1)
            plt.plot(1 - self.mu, 0, marker= '.', markersize= 15, color= self.b2['clr'], label= self.b2['name'])   # Plotting Smaller Body (Body 2)
            
            # Plotting Orbit of Spacecraft
            plt.plot( ys[0, 0], ys[0, 1], marker= '.', markersize= 10, color= 'magenta', label= f'{self.object_label}')
            plt.plot(ys[:, 0], ys[:, 1], color= 'orchid', label= f'Trajectory of {self.object_label}')

            plt.grid(linestyle= '--', linewidth= 0.2)
            plt.tight_layout()
            plt.title(self.title)
            plt.legend(bbox_to_anchor=(0.85, 0.8), loc='lower center')
            plt.show()
        
        elif self.plot == '3d':
            pass
            
        else:
            return ts, ys
