# ---------------------------- Groundtracks ------------------------------- #

import numpy as np
import matplotlib.pyplot as plt


city_list0 = [
	 'Seattle', 'Pasadena',                 # US
	 'New York', 'San Luis Obispo',
	 'Phoenix', 'Cape Canaveral',
	 'Mexico City', 'Villahermosa',         # Mexico
	 'New Delhi', 'Mumbai',                 # India
	 'Tirunelveli', 'Surat', 'Chennai',
	 'Olney', 'Norwich',                    # England
	 'Ponce',                               # Puerto Rico
	 'Berlin',                              # Germany
	 'Lyon',                                # France
	 'Vienna',                              # Austria
	 'Madrid', 'Sevilla', 'Barcelona',      # Spain
	 'Moscow',                              # Russia
	 'Baikonur',							# Kazakhstan
	 'Rome', 'Cortemaggiore',               # Italy
	 'Aalborg',                             # Denmark
	 'Sao Paulo',                           # Brazil
	 'Luxembourg City', 'Esch-sur-Alzette', # Luxembourg
	 'Toronto',                             # Canada
	 'Tokyo',                               # Japan
	 'Istanbul',                            # Turkey
	 'Jihlava',                             # Czech Republic
	 'Warsaw',                              # Poland
	 'Zagreb',                              # Croatia
	 'Sydney', 'Adelaide',                  # Australia
	 'Dubai',                               # UAE
	 'Port Louis',                          # Mauritius
	 'Casablanca',                          # Morocco
	 'Khartoum',                            # Sudan
	 'Tunis',                               # Tunisia
	 'Buenos Aires',                        # Argentina
	 'Cape Town',                           # South Africa
	 'Bucharest',                           # Romania
	 'Bogota',                              # Colombia
	 'Quito',                               # Ecuador
	 'Noordwijk',                           # Netherlands
	 'San Jose',                            # Costa Rica
	 'Stockholm',                           # Sweden
	 'Santiago',                            # Chile
	 'Jakarta',                             # Indonesia
	 'Antwerp',                             # Belgium
	 'Geneva',                              # Switzerland
	 'Manila',                              # Phillipines
	 'Porto', 'Ponta Delgada',              # Portugal
	 'Budapest',                            # Hungary
	 'Panama City',                         # Panama
	 'Cairo',                               # Egypt
	 'Seoul',                               # South Korea
	 'Broom Bridge',                        # Ireland
	 'Lima',                                # Peru
	 'Akure'                                # Nigeria
]

city_list1 = [
	 'Pasadena', 'Cape Canaveral',          # US
	 'Kodiak',
	 'Baikonur',							# Kazakhstan
	 'Kourou',								# French Guiana
	 
	 'New York', 'Houston',					# US
	 'Kansas City', 'Honolulu',
	 'Mexico City',         				# Mexico
	 'New Delhi', 'Mumbai',                 # India
	 'Ponce',                               # Puerto Rico
	 'Berlin',                              # Germany
	 'Madrid', 'Barcelona',      			# Spain
	 'Moscow',                              # Russia
	 'Sao Paulo',                           # Brazil
	 'Rome',               					# Italy
	 'Vancouver', 'Toronto',                # Canada
	 'Tokyo',                               # Japan
	 'Sydney', 'Adelaide',                  # Australia
	 'Dubai',                               # UAE
	 'Casablanca',                          # Morocco
	 'Khartoum',                            # Sudan
	 'Buenos Aires',                        # Argentina
	 'Cape Town',                           # South Africa
	 'Bogota',                              # Colombia
	 'San Jose',                            # Costa Rica
	 'Cairo',                               # Egypt
	 'Lagos',                               # Nigeria
	 'Seoul',                               # South Korea
	 'Bangkok',								# Thailand
	 'Singapore',							# Singapore
	 'Beijing',								# China
]

def city_dict():
    
	with open('world_cities.csv', 'r', errors= "ignore") as f:
		lines = f.readlines()

	header = lines[0]
	cities = {}

	for line in lines[1:]:
            
		line = line.split(',')

		# Creating a new dictionary for a given city
		try:
			cities[line[1]] = [float(line[2]), float(line[3])]
		except:
			pass

	return cities



def groundtracks(coordinates, city_names= city_list1, labels= None, cs= ['C9', 'c', 'C5', 'C1', 'C6', 'r', 'w', 'b', 'C6']):
    
	plt.style.use('dark_background')
	plt.figure()

	# Loading Coastline Coordinates
	coast_coords = np.genfromtxt('coastlines.csv', delimiter= ',')
	# Plotting Coastlines
	plt.plot(coast_coords[:, 0], coast_coords[:, 1], 'o', markersize= 0.3, color= 'magenta')

	# Plotting Orbits
	for n in range(len(coordinates)):
		if labels is None:
			label = str(n)
		else:
			label = labels[n]

		plt.plot(coordinates[n][0, 1], coordinates[n][0, 0], cs[n] + 'o', label= label)
		plt.plot(coordinates[n][1:, 1], coordinates[n][1:, 0], cs[n] + 'o', markersize= 0.75)

	# Plotting Cities
	cities = city_dict()
	n = 0

	for city in city_names:
		coord = cities[city]

		citycolor = 'white'
		if n <= 4:
			citycolor = 'lime'

		plt.plot([coord[1]], [coord[0]], 'o', markersize= 3, color= 'palegreen')

		# Alternating annotation above or below city
		if n % 2 == 0:
			xytext = (0, 3)
		else: 
			xytext = (0, -10)
		
		plt.annotate(city, [coord[1], coord[0]],
			   textcoords= 'offset points', xytext= xytext, ha= 'center', color= citycolor, fontsize= 'small')
		
		n += 1


	plt.grid(linestyle= 'dotted')
	plt.xlabel('Longitude (degrees $^\circ$)')
	plt.ylabel('Latitude (degrees $^\circ$)')
	plt.xlim(-180, 180)
	plt.ylim(-90, 90)
	plt.legend()
	plt.show()
