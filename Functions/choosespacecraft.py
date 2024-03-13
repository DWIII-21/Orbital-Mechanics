# ---------------- Function for Choosing Spacecraft --------------- #

#This function will call the chosen spacecraft's characteristics that are used in simulations

#OPTIONS FOR SPACECRAFT TYPE:
#Voyager 1(v)
#Orion(o)
#Starship(s)
#Falcon 9(f)


### INPUTS/OUTPUTS ###

# Inputs:
# Spacecraft Type

# Outputs:
# Name of Spacecraft (sc_type)
# Color that will be used to reference the spacecraft in a plot
# Color that will be used to reference the spacecraft's orbit in a plot
# DRY MASS of the spacecraft
# Specific Impulse of the spacecraft


###### CALCULATIONS ######


#Determing spacecraft type
def scchsr(sc):

    if sc == 'O' or sc == 'o' or sc == 'Orion':
        sc_type = 'Orion'
        scclr = 'skyblue'
        oclr = 'lightsteelblue'
        m_sc = 14045  #mass of Orion spacecraft in kg (from NASA; https://www.nasa.gov/pdf/166914main_FS_Orion508c.pdf)
        isp = 316 #Specific impulse in seconds (from Spaceflight101; https://spaceflight101.com/spacecraft/orion/)

    elif sc == 'S' or sc == 's' or sc == 'Starship':
        sc_type = 'Starship'
        scclr = 'mediumblue'
        oclr = 'silver'
        m_sc = 113398  #dry mass of SpaceX Starship spacecraft in kg (from Starship Weebly; https://spacex-guide.weebly.com/starship.html)
        isp = 327 #in seconds (from multiple sources)

    elif sc == 'V' or sc == 'v' or sc == 'Voyager':
        sc_type = 'Voyager 1'
        scclr = 'gold'
        oclr = 'gainsboro'
        m_sc = 721.9  #mass of Voyager spacecraft in kg (from NASA; https://nssdc.gsfc.nasa.gov/planetary/voyager.html)
        isp = 230 #in seconds (from NASA; https://ntrs.nasa.gov/api/citations/19660011758/downloads/19660011758.pdf)

    elif sc == 'F' or sc == 'f' or sc == 'Falcon 9' or sc == 'Falcon9':
        sc_type = 'Falcon 9'
        scclr = 'white'
        oclr = 'silver'
        m_sc = 549054  #mass of SpaceX Falcon 9 spacecraft in kg (from SpaceX; https://www.spacex.com/vehicles/falcon-9/)
        isp = 282 #in seconds (from Google)

    else:
        print("\nPlease restart and enter an available option for spacecraft type.")
        exit(1)

    return sc_type, scclr, oclr, m_sc, isp
