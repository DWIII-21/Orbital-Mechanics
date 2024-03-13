#-------------- Solar System Simulation ---------------#

import math as m
import pygame as pyg
import planorbchar as po
import numpy as np

pyg.init()

WIDTH, HEIGHT = 1000, 1000
WINDOW = pyg.display.set_mode((WIDTH,HEIGHT))
pyg.display.set_caption('Solar System Simulation')

FPS = 60

PLANET_FONT_SMALL = pyg.font.SysFont('calibri', 21)
PLANET_FONT_LARGE = pyg.font.SysFont('calibri', 30)
INSTRUCTIONS_FONT = pyg.font.SysFont('arial', 21)

[mu_S,eq_S,M_S,clr_S] = po.planchar('Sun')

class Instructions:

    def __init__(self,x,y,message):
        self.x = x
        self.y = y
        self.message = message

        self.instructions_text = INSTRUCTIONS_FONT.render(self.message, 1, 'White')

    def draw(self):
        WINDOW.blit(self.instructions_text, (self.x, self.y))

class Body:
    AU = 1.496*10**11  # 1 AU in m
    G = 6.67428*10**-11  # Gravitational Constant
    SCALE = 300/AU  # 1 AU = 100 pixels
    TIMESTEP = 3600*24 # 1 Day
    line_thickness = 2  # in number of pixels

    def __init__(self, au, radius, name): # Where au is the number of AU that the planet is from the Sun
        self.x = au * self.AU
        self.y = 0
        self.radius = radius
        self.name = name
        [self.mu,self.eq,self.mass,self.color] = po.planchar(name)
        [self.temp,self.moons,self.gravity,self.inclination,self.day_length] = po.planmisc(name)
        if not name == 'Sun':
            [self.a,self.e,self.rp,self.ra,self.b,self.p,self.h] = po.orbchar(name)

        self.orbit = []
        self.sun = False  # Tells us if the body in question is the Sun
        self.distance_to_sun = 0

        self.x_vel = 0
        self.y_vel = 0

        # Finding average velocity of the planet in m/s
        if self.x != 0:
            self.avg_vel = m.sqrt(self.G*M_S/abs(self.x))
        if self.x < 0:
            self.y_vel = self.avg_vel
        elif self.x > 0:
            self.y_vel = -self.avg_vel

        # Year Length Calculation
        self.year_length = m.sqrt(abs(self.x)**3 * (4*m.pi**2)/(self.G*(self.mass + M_S)))/(365*86477)

    def attraction(self,body2):
        force_x = 0
        force_y = 0
        if not self.sun:
            body2_x, body2_y = body2.x, body2.y
            distance_x = body2_x - self.x
            distance_y = body2_y - self.y
            distance = m.sqrt(distance_x**2 + distance_y**2)
            self.distance = distance

            self.radii = []
            self.radii.append(self.distance)

            if body2.sun:  # Check if other body in question is the Sun
                self.distance_to_sun = distance

            force = self.G*self.mass*body2.mass/distance**2
            theta = m.atan2(distance_y,distance_x)
            force_x = m.cos(theta)*force
            force_y = m.sin(theta)*force
        
        return force_x, force_y
    
    def update_position(self,bodies):
        total_Fx = total_Fy = 0
        
        for body in bodies:
            if self == body:
                continue

            Fx, Fy = self.attraction(body)
            total_Fx += Fx
            total_Fy += Fy

        self.x_vel += total_Fx / self.mass * self.TIMESTEP  # From F = m*a
        self.y_vel += total_Fy / self.mass * self.TIMESTEP

        self.x += self.x_vel * self.TIMESTEP
        self.y += self.y_vel * self.TIMESTEP
        self.orbit.append((self.x,self.y))
    
    def draw_bodies(self,window,display_message,numberpressed_message):
        x = self.x * self.SCALE + WIDTH/2
        y = self.y * self.SCALE + HEIGHT/2
 
        if len(self.orbit) > 2 and self.sun == False:
            updated_points = []
            for point in self.orbit:
                    x, y = point
                    x = x * self.SCALE + WIDTH/2
                    y = y * self.SCALE + HEIGHT/2
                    updated_points.append((x,y))

            pyg.draw.lines(window, self.color, False, updated_points, self.line_thickness) # Line thickness of 2 pixels

        pyg.draw.circle(window,self.color,(x,y),self.radius)

        text = PLANET_FONT_LARGE.render(display_message, 1 , 'lime')
        window.blit(text, ((x - text.get_width()/2, y - text.get_height()/2 + 20)))
        numberpressed_text = PLANET_FONT_LARGE.render(numberpressed_message, 1, 'Green')
        window.blit(numberpressed_text,((21, HEIGHT - 75) ))

# Instructions
cycle = Instructions(5, 5, "To Cycle through Planetary Stats:")
hold = Instructions(15, 30, "HOLD ANY NUMBER KEY (0-9)")
inout = Instructions(WIDTH/2 + 70, 5, "To Cycle Between Inner and Outer Planets:")
space = Instructions(WIDTH - 300, 30, "PRESS SPACE")

def main():
    clock = pyg.time.Clock()

    # Creating the Sun
    sun = Body(0, 30, 'Sun')
    sun.sun = True
    # Creating Mercury
    mercury = Body(0.387098, 8, 'Mercury')
    # Creating Venus
    venus = Body(0.723328, 14, 'Venus')
    # Creating Earth
    earth = Body(-1, 16, 'Earth')
    # Creating Mars
    mars = Body(-1.523706, 12, 'Mars')
    # Creating Jupiter
    jupiter = Body(5.20149, 18, 'Jupiter')
    # Creating Saturn
    saturn = Body(9.54327, 16, 'Saturn')
    # Creating Uranus
    uranus = Body(-19.17113, 11, 'Uranus')
    # Creating Neptune
    neptune = Body(-29.99375, 10, 'Neptune')
    # Creating Pluto
    pluto = Body(39.2305, 3, 'Pluto')

    bodies = [sun,mercury,venus,earth,mars]

    outer = False
    run = True
    while run:
        clock.tick(FPS)
        WINDOW.fill('Black')
        number_of_message = 0
        display_message = ''
        numberpressed_message = 'N/A'

        cycle.draw()
        hold.draw()
        inout.draw()
        space.draw()
        
        keys = pyg.key.get_pressed()
        for body in bodies:
            if keys[pyg.K_1] or keys[pyg.K_KP_1]:
                number_of_message = 1
            elif keys[pyg.K_2] or keys[pyg.K_KP_2]:
                number_of_message = 2
            elif keys[pyg.K_3] or keys[pyg.K_KP_3]:
                number_of_message = 3
            elif keys[pyg.K_4] or keys[pyg.K_KP_4]:
                number_of_message = 4
            elif keys[pyg.K_5] or keys[pyg.K_KP_5]:
                number_of_message = 5
            elif keys[pyg.K_6] or keys[pyg.K_KP_6]:
                number_of_message = 6
            elif keys[pyg.K_7] or keys[pyg.K_KP_7]:
                number_of_message = 7
            elif keys[pyg.K_8] or keys[pyg.K_KP_8]:
                number_of_message = 8
            elif keys[pyg.K_9] or keys[pyg.K_KP_9]:
                number_of_message = 9
            elif keys[pyg.K_0] or keys[pyg.K_KP_0]:
                number_of_message = 10
                

        if keys[pyg.K_SPACE]:
            if outer == False:
                sun.x = 0
                sun.y = 0
                Body.SCALE = 12/Body.AU
                Body.line_thickness = 1
                bodies = [sun,jupiter,saturn,uranus,neptune,pluto]
                for i in range(len(bodies)):
                    if bodies[i].sun == False:
                        Body.TIMESTEP = 3600*24*50
                outer = True
                clock.tick(30)
            elif outer == True:
                sun.x = 0
                sun.y = 0
                Body.SCALE = 300/Body.AU
                Body.line_thickness = 2
                for i in range(len(bodies)):
                    if bodies[i].sun == False:
                        Body.TIMESTEP = 3600*24
                bodies = [sun,mercury,venus,earth,mars]
                outer = False
                clock.tick(FPS)

        for event in pyg.event.get():
            if event.type == pyg.QUIT:
                run = False
                break

        for body in bodies:
            if number_of_message == 1:
                display_message = f'{body.name}'
                numberpressed_message = 'Planet Names'
            elif number_of_message == 2:
                if not body.sun:
                    display_message = f'{round(body.distance_to_sun/1000,1)} km'
                numberpressed_message = 'Distance to Sun'
            elif number_of_message == 3:
                if not body.sun:
                    display_message = f'≈ {round(body.year_length,2)} Years'
                numberpressed_message = 'Orbital Period'
            elif number_of_message == 4:
                if not body.sun:
                    display_message = f'≈ {round(body.e,3)}'
                numberpressed_message = 'Orbit Eccentricity'
            elif number_of_message == 5:
                display_message = f'{body.temp}°F'
                numberpressed_message = 'Avg. Temperature'
            elif number_of_message == 6:
                if not body.sun:
                    if body.moons == 1:
                        display_message = f'{body.moons} Moon'
                    else:
                        display_message = f'{body.moons} Moons'
                numberpressed_message = 'Number of Moons'  
            elif number_of_message == 7:
                display_message = f'{round(body.eq*2,2)} km'
                numberpressed_message = 'Planet Diameter'
            elif number_of_message == 8:
                if not body.sun:
                    display_message = f'{body.gravity} m/s^2'
                numberpressed_message = 'Planet Gravity'
            elif number_of_message == 9:
                if not body.sun:
                    display_message = f'{body.inclination}°'
                numberpressed_message = 'Orbit Inclination'
            elif number_of_message == 10:
                if not body.sun:
                    display_message = f'{body.day_length} Hours'
                numberpressed_message = 'Planet Day Length'

            body.update_position(bodies)
            body.draw_bodies(WINDOW,display_message,numberpressed_message)
        

        # Displaying label for which planet set is being displayed
        planetdisplay_message = 'Currently Displaying:'
        planetdisplay_text = INSTRUCTIONS_FONT.render(planetdisplay_message,1,'White')
        if outer:
            planetset_message = 'OUTER PLANETS'
            planetset_text = PLANET_FONT_LARGE.render(planetset_message, 1, 'Orange')
        else:
            planetset_message = 'INNER PLANETS'
            planetset_text = PLANET_FONT_LARGE.render(planetset_message,1,'Yellow')
        WINDOW.blit(planetdisplay_text, (21, HEIGHT - 100))
        WINDOW.blit(planetset_text, (WIDTH - 250, HEIGHT - 75))
        
        pyg.display.update()
    
    pyg.quit()

if __name__ == '__main__':
    main()
