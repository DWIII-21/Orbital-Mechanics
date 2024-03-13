#-----------  Earth to Mars Heliocentric Transfer -----------#

import math as m
import matplotlib.pyplot as plt
import numpy as np
import Julian as jul
import Gregorian as g
import ephem as eph
import class2cart as cl
import class2cartcomponents as cm
import lambert as lmb
import C3calc as c3
import f_and_g as fg
import planorbchar as po

def main():

    # Establishing initial variables
    muSun = 1.327*10**11
    mu = muSun
    AU_2_km = 1.469*10**8  # AU to km conversion
    aE = 1  # Semi-major axis of Earth's orbit in AU
    aM = 1.523706  # Semi-major axis of Mars' orbit in AU
    e_Earth = 0.016882 # Earth's Orbit's Eccentricity
    e_Mars = 0.09352 # Mars' Orbit's Eccentricity
    size = 10
    h = 100  # Altitude in km above Mars at which entry velocity is being measured

    [mu_Earth, eq_Earth, M_Earth, clr_Earth] = po.planchar('Earth')
    [mu_Mars, eq_Mars, M_Mars, clr_Mars] = po.planchar('Mars')

    # Calculating Synodic Period
    n1 = (mu/((aE*AU_2_km)**3))**(1/2)
    n2 = (mu/((aM*AU_2_km)**3))**(1/2)
    tsynd = (2*m.pi)/(n1-n2)/(60*60*24*365)  # Synodic period in years

    # Setting a minimum and maximum time of flight (in days)
    TOFmin = 30
    TOFmax = TOFmin + tsynd*365
    TOFvec = np.linspace(TOFmin,TOFmax,size)

    # Establish range for minimum and maximum Departure and Arrival Dates
    JDdpmin = jul.julianpre(2023,11,21,0,0)
    JDdpmax = JDdpmin + tsynd*365

    # Establishing meshgrids for departure and arrival
    JDdpvec = np.linspace(JDdpmin,JDdpmax,size)
    [dpmesh,TOF] = np.meshgrid(JDdpvec,TOFvec)
    armesh = dpmesh + TOF
    
    # Initializing C3 and v_entry matrices and setting maximum values for each measurement
    C3 = np.zeros((len(JDdpvec),len(TOFvec)),dtype=float)
    v_entry = np.zeros((len(JDdpvec),len(TOFvec)),dtype=float)
    C3_max = 500
    v_entry_max = 50

    # Indexing values into C3 and v_entry matrices
    for i in range(len(JDdpvec)):

        print(f'Loading C3 and Entry Velocity Matrices. Status: {int(int(i)/int(len(JDdpvec))*100)} %')  # Display a status until code runs to completion

        for j in range(len(JDdpvec)):
            
            [C3_check, v_entry_check] = c3.C3calc(dpmesh[i,j],armesh[i,j],mu,'Earth','Mars',h)

            # Eliminating outlier values that exceed general range of calculations for C3 and v_entry
            if C3_check >= C3_max or v_entry_check >= v_entry_max:
                if v_entry_check < v_entry_max and C3_check >= C3_max:
                    C3[i,j] = C3_max
                    v_entry[i,j] = v_entry_check
                elif C3_check < C3_max and v_entry_check >= v_entry_max:
                    C3[i,j] = C3_check
                    v_entry[i,j] = v_entry_max
                else:
                    C3[i,j] = C3_max
                    v_entry[i,j] = v_entry_max

            else:    
                [C3[i,j],v_entry[i,j]] = c3.C3calc(dpmesh[i,j],armesh[i,j],mu,'Earth','Mars',h)
    
    print('Loading C3 and Entry Velocity Matrices. Status: 100 %')
    print('\nDone Loading C3 and Entry Velocity Matrices')
    print('Completed Step 1/4\n')

    C3min = np.min(C3)
    C3index_min = [0,0]
    for i in range(len(JDdpvec)):
        print(f'Finding Minimum C3. Status: {int(int(i)/int(len(JDdpvec))*100)} %')  # Display a status until code runs to completion
        for j in range(len(JDdpvec)):
            # Checking for location of minimum values in each matrix
            if C3[i,j] == C3min:
                C3index_min = [i,j]
    
    print('Finding Minimum C3. Status: 100 %')
    print('\nDone finding location of minimum C3')
    print('Completed Step 2/4\n')        

    # Finding Departure Date and Time of Flight with lowest C3
    DP_C3min = dpmesh[C3index_min[0],C3index_min[1]]
    [DP_C3min_month, DP_C3min_day, DP_C3min_year] = g.gregorian(DP_C3min)
    AR_C3min = armesh[C3index_min[0],C3index_min[1]]
    [AR_C3min_month, AR_C3min_day, AR_C3min_year] = g.gregorian(AR_C3min)
    Chosen_TOF_Days = AR_C3min - DP_C3min 

    launch_window_min_DP, launch_window_max_DP = DP_C3min - 15, DP_C3min + 15
    [month_min_DP,day_min_DP,year_min_DP],[month_max_DP,day_max_DP,year_max_DP] = g.gregorian(launch_window_min_DP), g.gregorian(launch_window_max_DP)

    launch_window_min_AR, launch_window_max_AR = 2.4609*10**6, 2.461025*10**6  # Values chosen from observation of the pork chop plot
    window_dpvec = np.linspace(launch_window_min_DP,launch_window_max_DP,size)
    window_arvec = np.linspace(launch_window_min_AR,launch_window_max_AR,size)
    window_v_entry = np.zeros((len(window_dpvec),len(window_arvec)),dtype=float)

    blank = 0
    for i in range(len(window_dpvec)):

        print(f'Loading Launch Window Matix. Status: {int(int(i)/int(len(window_dpvec))*100)} %')  # Display a status until code runs to completion

        for j in range(len(window_arvec)):
            
            blank, window_v_entry[i,j] = c3.C3calc(window_dpvec[i],window_arvec[j],mu,'Earth','Mars',h)

    print('Loading Launch Window Matix. Status: 100 %')
    print('\nDone Loading Launch Window Matrix')
    print('Completed Step 3/4\n')
    
    v_entrymin = np.min(window_v_entry)
    v_entrymax = np.max(window_v_entry)
    v_entryindex_min = [0.0]
    v_entryindex_max = [0,0]
    for i in range(len(window_dpvec)):
        print(f'Finding Minimum and Maximum Entry Velocities. Status: {int(int(i)/int(len(window_dpvec))*100)} %')  # Display a status until code runs to completion
        for j in range(len(window_dpvec)):
            # Checking for location of minimum and maximum entry velocity in the matrix
            if window_v_entry[i,j] == v_entrymin:
                v_entryindex_min = [i,j]
            if window_v_entry[i,j] == v_entrymax:
                v_entryindex_max = [i,j]

    print('Finding Minimum and Maximum Entry Velocities. Status: 100 %')
    print('Completed Step 4/4\n')

    # Data required to plot trajectory to Mars

    Chosen_JDdp = DP_C3min
    Chosen_JDar = AR_C3min

    EarthDP_a, EarthDP_e, EarthDP_i, EarthDP_omega, EarthDP_w, EarthDP_M = eph.ephem(Chosen_JDdp,'Earth')
    EarthAR_a, EarthAR_e, EarthAR_i, EarthAR_omega, EarthAR_w, EarthAR_M = eph.ephem(Chosen_JDar,'Earth')
    MarsDP_a, MarsDP_e, MarsDP_i, MarsDP_omega, MarsDP_w, MarsDP_M = eph.ephem(Chosen_JDdp,'Mars')
    MarsAR_a, MarsAR_e, MarsAR_i, MarsAR_omega, MarsAR_w, MarsAR_M = eph.ephem(Chosen_JDar,'Mars')

    # Finding each planet's classical coordinates at depature and arrival (captial R and V signify vector measurements of these values)
    R_EarthDP, V_EarthDP = cl.class2cart(mu, EarthDP_a, EarthDP_e, EarthDP_i, EarthDP_omega, EarthDP_w, EarthDP_M)
    R_EarthAR, V_EarthAR = cl.class2cart(mu, EarthAR_a, EarthAR_e, EarthAR_i, EarthAR_omega, EarthAR_w, EarthAR_M)
    R_MarsDP, V_MarsDP = cl.class2cart(mu, MarsDP_a, MarsDP_e, MarsDP_i, MarsDP_omega, MarsDP_w, MarsDP_M)
    R_MarsAR, V_MarsAR = cl.class2cart(mu, MarsAR_a, MarsAR_e, MarsAR_i, MarsAR_omega, MarsAR_w, MarsAR_M)


    # Rearranging velocity vectors to allow proper vector math, cannot change radius because of its use in the lambert function
    V_EarthDP = V_EarthDP.reshape((1,3))
    V_EarthAR = V_EarthAR.reshape((1,3))
    V_MarsDP = V_MarsDP.reshape((1,3))
    V_MarsAR = V_MarsAR.reshape((1,3))

    Chosen_TOF = (Chosen_JDar - Chosen_JDdp)*24*3600  # TOF in seconds

    # Finding semi-major axis and parameter of transfer orbit
    [a_trsf, p_trsf, VDP_trsf,VAR_trsf,conv] = lmb.lambert(R_EarthDP, R_MarsAR, Chosen_TOF, mu, 1, 100, 10**(-14), 100)

    # Finding scalar magnitudes of each vector (denoted by lowercase letters)
    r_EarthDP = np.linalg.norm(R_EarthDP)
    r_EarthAR = np.linalg.norm(R_EarthAR)
    r_MarsDP = np.linalg.norm(R_MarsDP)
    r_MarsAR = np.linalg.norm(R_MarsAR)

    v_EarthDP = np.linalg.norm(V_EarthDP)
    v_EarthAR = np.linalg.norm(V_EarthAR)
    v_MarsDP = np.linalg.norm(V_MarsDP)
    v_MarsAR = np.linalg.norm(V_MarsAR)

    e_trsf = m.sqrt(1-p_trsf/a_trsf) # Eccentricity of the transfer orbit
    p_Earth = aE * AU_2_km  # Parameter of Earth's orbit about the Sun
    p_Mars = aM * AU_2_km  # Parameter of Mars' orbit about the Sun

    # Finding true anomaly values at depature and arrival for each orbit
    nu_EarthDP = m.acos((p_Earth/r_EarthDP-1) * 1/e_Earth)
    nu_EarthAR = m.acos((p_Earth/r_EarthAR-1) * 1/e_Earth)
    nu_MarsDP = m.acos((p_Mars/r_MarsDP-1) * 1/e_Mars)
    nu_MarsAR = m.acos((p_Mars/r_MarsAR-1) * 1/e_Mars)
    nu_trsfDP = m.acos((p_trsf/r_EarthDP-1) * 1/e_trsf)
    nu_trsfAR = m.acos((p_trsf/r_MarsAR-1) * 1/e_trsf)

    fg_Earth = fg.fg(0, 2*np.pi, e_Earth, R_EarthDP, p_Earth, V_EarthDP, {'size': 5*10**4})
    fg_Mars = fg.fg(0, 2*np.pi, e_Mars, R_MarsDP, p_Mars, V_MarsDP, {'size': 5*10**4})
    fg_trsf = fg.fg(0, 2*np.pi, e_trsf,  R_EarthDP, p_trsf, VDP_trsf, {'size': 5*10**4})

    # Printing results for ideal launch window
    print(f"\nThe mission yields a minimum C3 value of {round(C3min,2)} km^2/s^2 for the given minimum launch date.")
    print(f"\nThis places the ideal 30-day launch window between:\n"
          f"{month_min_DP}/{day_min_DP}/{year_min_DP} and {month_max_DP}/{day_max_DP}/{year_max_DP}")
    print(f"\nFor this launch window, the minimum and maximum entry velocities into Mars at an altitude of {h} km would be:\n"
          f"{round(v_entrymin,2)} km/s and {round(v_entrymax,2)} km/s, respectively.\n")
    # Printing Chosen Departure and Time of Flight
    print(f"If the spacecraft departs on the date which requires the lowest C3,\nit would leave on {DP_C3min_month}/{DP_C3min_day}/{DP_C3min_year} and arrive on {AR_C3min_month}/{AR_C3min_day}/{AR_C3min_year},\ngiving the mission a time of flight of about {int(Chosen_TOF_Days)} days.\n")

    ##################
    # UNEDITED PLOTS #
    ##################

    # Graphing Contour Plots for C3 and vinfm

    # Level for each contour plot
    C3_levels = np.linspace(20,500,20)
    v_entry_levels = np.linspace(1,50)
    lw = 1 # Setting linewidth for launch window

    # Setting min and max values for axes of secondary plots
    xmin = launch_window_min_DP - 10
    xmax = launch_window_max_DP + 10
    ymin = launch_window_min_AR - 10
    ymax = launch_window_max_AR + 10

    # Calculating values necessary for launch window
    Absolute_AR_min = float(np.min(armesh))
    Absolute_AR_max = float(np.max(armesh))
    plot_height = Absolute_AR_max - Absolute_AR_min
    AR_min_height = launch_window_min_AR - Absolute_AR_min
    AR_max_height = launch_window_max_AR - Absolute_AR_min
    AR_min_ratio = AR_min_height/plot_height
    AR_max_ratio = AR_max_height/plot_height

    # Calculating values necessary for launch window in Zoomed Plot
    zoomed_height = ymax - ymin
    zoomed_min_height = launch_window_min_AR - ymin
    zoomed_max_height = launch_window_max_AR - ymin
    zoomed_min_ratio = zoomed_min_height/zoomed_height
    zoomed_max_ratio = zoomed_max_height/zoomed_height

    plt.subplot(121)
    plt.title('Characteristic Energy (C3) Leaving Earth')
    plt.xlabel('Julian Departure Dates')
    plt.ylabel('Julian Arrival Dates')
    C3plot = plt.contour(dpmesh,armesh,C3,levels=C3_levels)
    plt.clabel(C3plot, inline=1, fontsize=8.5)
    cbar1 = plt.colorbar(C3plot)
    cbar1.ax.set_ylabel('C3 (in km^2/s^2)')
    plt.grid(color = 'lightgray', linestyle='--',linewidth = 0.3)

    plt.subplot(122)
    plt.title(f"Entry Speed into Mars' Atmosphere\nat Altitude of {h} km")
    plt.xlabel('Julian Departure Dates')
    plt.ylabel('Julian Arrival Dates')
    v_entryplot = plt.contour(dpmesh,armesh,v_entry,levels=v_entry_levels)
    plt.clabel(v_entryplot, inline=1, fontsize=9)
    cbar2 = plt.colorbar(v_entryplot)
    cbar2.ax.set_ylabel('Entry Speed (in km/s)')
    plt.grid(color = 'lightgray', linestyle='--',linewidth = 0.3)

    plt.show()

    ############
    # C3 PLOTS #
    ############

    plt.subplot(121)
    plt.title('Characteristic Energy (C3) Leaving Earth')
    plt.xlabel('Julian Departure Dates')
    plt.ylabel('Julian Arrival Dates')
    C3plot = plt.contour(dpmesh,armesh,C3,levels=C3_levels)
    plt.clabel(C3plot, inline=1, fontsize=8.5)
    cbar1 = plt.colorbar(C3plot)
    cbar1.ax.set_ylabel('C3 (in km^2/s^2)')
    plt.axvline(launch_window_min_DP,color='red',linestyle='--',linewidth = lw)
    plt.axvline(launch_window_max_DP,color='red',linestyle='--',linewidth = lw)
    plt.axhline(launch_window_min_AR,color='red',linestyle='--',linewidth = lw)
    plt.axhline(launch_window_max_AR,color='red',linestyle='--',linewidth = lw)
    plt.axvspan(launch_window_min_DP, launch_window_max_DP, AR_min_ratio, AR_max_ratio, color='red',alpha=0.5,label='Launch Window')

    plt.legend(bbox_to_anchor=(0.15, 0.9), loc='lower center')
    plt.grid(color = 'lightgray', linestyle='--',linewidth = 0.3)

    ax1 = plt.subplot(122)
    plt.title('Characteristic Energy (C3) Leaving Earth (ZOOMED IN)')
    plt.xlabel('Julian Departure Dates')
    plt.ylabel('Julian Arrival Dates')
    C3plot = plt.contour(dpmesh,armesh,C3,levels=C3_levels)
    plt.clabel(C3plot, inline=1, fontsize=8.5)
    ax1.set_xlim((xmin,xmax))
    ax1.set_ylim((ymin,ymax))
    cbar1 = plt.colorbar(C3plot)
    cbar1.ax.set_ylabel('C3 (in km^2/s^2)')
    # Plotting location of min C3
    plt.plot(dpmesh[C3index_min[0],C3index_min[1]],armesh[C3index_min[0],C3index_min[1]],marker='o',markersize=7,color='lime',label='Location of\nMinimum C3')
    # Plotting launch window
    plt.axvline(launch_window_min_DP,color='red',linestyle='--',linewidth = lw)
    plt.axvline(launch_window_max_DP,color='red',linestyle='--',linewidth = lw)
    plt.axhline(launch_window_min_AR,color='red',linestyle='--',linewidth = lw)
    plt.axhline(launch_window_max_AR,color='red',linestyle='--',linewidth = lw)
    plt.axvspan(launch_window_min_DP, launch_window_max_DP, zoomed_min_ratio, zoomed_max_ratio, color='red',alpha=0.25,label='Launch Window')
    #plt.axvspan(launch_window_min_DP, launch_window_max_DP, 0.07, 0.93, color='red',alpha=0.25,label='Launch Window')

    plt.legend(bbox_to_anchor=(0.17, 0.89), loc='lower center')
    plt.grid(color = 'lightgray', linestyle='--',linewidth = 0.3)

    plt.show()

    ########################
    # ENTRY VELOCITY PLOTS #
    ########################

    plt.subplot(121)
    plt.title(f"Entry Speed into Mars' Atmosphere\nat Altitude of {h} km")
    plt.xlabel('Julian Departure Dates')
    plt.ylabel('Julian Arrival Dates')
    v_entryplot = plt.contour(dpmesh,armesh,v_entry,levels=v_entry_levels)
    plt.clabel(v_entryplot, inline=1, fontsize=9)
    cbar2 = plt.colorbar(v_entryplot)
    cbar2.ax.set_ylabel('Entry Speed (in km/s)')
    # Plotting launch window
    plt.axvline(launch_window_min_DP,color='red',linestyle='--',linewidth = lw)
    plt.axvline(launch_window_max_DP,color='red',linestyle='--',linewidth = lw)
    plt.axhline(launch_window_min_AR,color='red',linestyle='--',linewidth = lw)
    plt.axhline(launch_window_max_AR,color='red',linestyle='--',linewidth = lw)
    plt.axvspan(launch_window_min_DP, launch_window_max_DP, AR_min_ratio, AR_max_ratio, color='red',alpha=0.5,label='Launch Window')

    plt.legend(bbox_to_anchor=(0.15, 0.85), loc='lower center')
    plt.grid(color = 'lightgray', linestyle='--',linewidth = 0.3)

    ax2 = plt.subplot(122)
    plt.title(f"Entry Speed into Mars' Atmosphere\nat Altitude of {h} km (ZOOMED IN)")
    plt.xlabel('Julian Departure Dates')
    plt.ylabel('Julian Arrival Dates')
    v_entryplot = plt.contour(dpmesh,armesh,v_entry,levels=v_entry_levels)
    plt.clabel(v_entryplot, inline=1, fontsize=9)
    ax2.set_xlim((xmin,xmax))
    ax2.set_ylim((ymin,ymax))
    cbar2 = plt.colorbar(v_entryplot)
    cbar2.ax.set_ylabel('Entry Speed (in km/s)')
    # Plotting location of min and max entry velo
    plt.plot(window_dpvec[v_entryindex_min[0]],window_arvec[v_entryindex_min[1]],marker='o',markersize=7,color='lime',label='Location of\nMinimum Entry Velocity')
    plt.plot(window_dpvec[v_entryindex_max[0]],window_arvec[v_entryindex_max[1]],marker='o',markersize=7,color='red',label='Location of\nMaximum Entry Velocity')
    # Plotting launch window
    plt.axvline(launch_window_min_DP,color='red',linestyle='--',linewidth = lw)
    plt.axvline(launch_window_max_DP,color='red',linestyle='--',linewidth = lw)
    plt.axhline(launch_window_min_AR,color='red',linestyle='--',linewidth = lw)
    plt.axhline(launch_window_max_AR,color='red',linestyle='--',linewidth = lw)
    plt.axvspan(launch_window_min_DP, launch_window_max_DP, zoomed_min_ratio, zoomed_max_ratio, color='red',alpha=0.25,label='Launch Window')

    plt.legend(bbox_to_anchor=(0.21, 0.85), loc='lower center')
    plt.grid(color = 'lightgray', linestyle='--',linewidth = 0.3)

    plt.show()

    ##################################
    # PROPAGATING TRAJECTORY TO MARS #
    ##################################

    '''
    This section propagates a spacecraft's trajectory from Earth to Mars on a departure date of 10/23/24. 
    This date is found with in the sweet spot calculated above of 10/8/24 and 11/7/24. 
    '''

    Chosen_JDdp = jul.julianpre(2024, 10, 23, 0, 0)
    Chosen_JDar = Chosen_JDdp + 367

    EarthDP_a, EarthDP_e, EarthDP_i, EarthDP_omega, EarthDP_w, EarthDP_M = eph.ephem(Chosen_JDdp,'Earth')
    EarthAR_a, EarthAR_e, EarthAR_i, EarthAR_omega, EarthAR_w, EarthAR_M = eph.ephem(Chosen_JDar,'Earth')
    MarsDP_a, MarsDP_e, MarsDP_i, MarsDP_omega, MarsDP_w, MarsDP_M = eph.ephem(Chosen_JDdp,'Mars')
    MarsAR_a, MarsAR_e, MarsAR_i, MarsAR_omega, MarsAR_w, MarsAR_M = eph.ephem(Chosen_JDar,'Mars')

    # Finding each planet's classical coordinates at depature and arrival (captial R and V signify vector measurements of these values)
    R_EarthDP, V_EarthDP = cl.class2cart(mu, EarthDP_a, EarthDP_e, EarthDP_i, EarthDP_omega, EarthDP_w, EarthDP_M)
    R_EarthAR, V_EarthAR = cl.class2cart(mu, EarthAR_a, EarthAR_e, EarthAR_i, EarthAR_omega, EarthAR_w, EarthAR_M)
    R_MarsDP, V_MarsDP = cl.class2cart(mu, MarsDP_a, MarsDP_e, MarsDP_i, MarsDP_omega, MarsDP_w, MarsDP_M)
    R_MarsAR, V_MarsAR = cl.class2cart(mu, MarsAR_a, MarsAR_e, MarsAR_i, MarsAR_omega, MarsAR_w, MarsAR_M)

    # Rearranging velocity vectors to allow proper vector math, cannot change radius because of its use in the lambert function
    V_EarthDP = V_EarthDP.reshape((1,3))
    V_EarthAR = V_EarthAR.reshape((1,3))
    V_MarsDP = V_MarsDP.reshape((1,3))
    V_MarsAR = V_MarsAR.reshape((1,3)) 

    Chosen_TOF = (Chosen_JDar - Chosen_JDdp)*24*3600  # TOF in seconds

    # Finding semi-major axis and parameter of transfer orbit
    [a_trsf, p_trsf, VDP_trsf,VAR_trsf,conv] = lmb.lambert(R_EarthDP, R_MarsAR, Chosen_TOF, mu, 1, 100, 10**(-14), 100)

    # Finding scalar magnitudes of each vector (denoted by lowercase letters)
    r_EarthDP = np.linalg.norm(R_EarthDP)
    r_EarthAR = np.linalg.norm(R_EarthAR)
    r_MarsDP = np.linalg.norm(R_MarsDP)
    r_MarsAR = np.linalg.norm(R_MarsAR)

    e_trsf = m.sqrt(1 - p_trsf/a_trsf) # Eccentricity of the transfer orbit
    p_Earth = aE * AU_2_km * (1 - e_Earth**2)  # Parameter of Earth's orbit about the Sun
    p_Mars = aM * AU_2_km * (1 - e_Mars**2)  # Parameter of Mars' orbit about the Sun

    # Finding true anomaly values at depature and arrival for each orbit
    nu_EarthDP = m.acos((p_Earth/r_EarthDP-1) * 1/e_Earth)
    nu_EarthAR = m.acos((p_Earth/r_EarthAR-1) * 1/e_Earth)
    nu_MarsDP = m.acos((p_Mars/r_MarsDP-1) * 1/e_Mars)
    nu_MarsAR = m.acos((p_Mars/r_MarsAR-1) * 1/e_Mars)

    nu_Earth = np.linspace(0, 2*np.pi, 100)
    nu_Mars = np.linspace(0, 2*np.pi, 100)

    R_Earth = np.zeros((100,3))
    R_Mars = np.zeros((100,3))

    R_EarthDP = R_EarthDP.reshape((1,3))
    R_MarsDP = R_MarsDP.reshape((1,3))
    for b in range(100):

        f_Earth = fg.lamF(p_Earth, e_Earth, nu_Earth[b], nu_EarthDP)
        f_Mars = fg.lamF(p_Mars, e_Mars, nu_Mars[b], nu_MarsDP)
        g_Earth = fg.lamG(p_Earth, e_Earth, nu_Earth[b], nu_EarthDP)
        g_Mars = fg.lamG(p_Mars, e_Mars, nu_Mars[b], nu_MarsDP)

        R_Earth[b,0], R_Earth[b,1], R_Earth[b,2] = (f_Earth*R_EarthDP + g_Earth*V_EarthDP)[0,0], (f_Earth*R_EarthDP + g_Earth*V_EarthDP)[0,1], (f_Earth*R_EarthDP + g_Earth*V_EarthDP)[0,2]
        R_Mars[b,0], R_Mars[b,1], R_Mars[b,2] = (f_Mars*R_MarsDP + g_Mars*V_MarsDP)[0,0], (f_Mars*R_MarsDP + g_Mars*V_MarsDP)[0,1], (f_Mars*R_MarsDP + g_Mars*V_MarsDP)[0,2]

    R_EarthDP = R_EarthDP.reshape((3,1))
    R_MarsDP = R_MarsDP.reshape((3,1))

    # Propagating Earth and Mars Trajectories (Using F & G Functions)
    fg_Earth = fg.fg(nu_EarthDP, nu_EarthDP + 2*np.pi, e_Earth, R_EarthDP, p_Earth, V_EarthDP, {'size': 5*10**4})
    fg_Mars = fg.fg(nu_MarsAR, nu_MarsAR + 2*np.pi, e_Mars, R_MarsAR, p_Mars, V_MarsAR, {'size': 5*10**4})

    # Propagating Transfer Trajectory (Using RK4 Solver)
    y0= [R_EarthDP[0][0], R_EarthDP[1][0], R_EarthDP[2][0], VDP_trsf[0], VDP_trsf[1], VDP_trsf[2]]  # Initial position of S/C on transfer orbit
    transfer_orbit = prop.rk4_2Body(cb= 'sun', dt= 5*10**5, tspan= 3.15*10**7, show= False, \
                                        y0= y0)

    R_EarthDP = R_EarthDP.reshape((3,1))
    R_MarsDP = R_MarsDP.reshape((3,1))

    ###############################
    # PLOTTING TRAJECTORY TO MARS #
    ###############################

    plt.style.use('dark_background')
    ax = plt.axes(projection = '3d')
    ax.plot(fg_Earth[:,0],fg_Earth[:,1],fg_Earth[:,2], color = 'blue', label = "Earth's Orbit")
    ax.plot(fg_Mars[:,0],fg_Mars[:,1],fg_Mars[:,2], color = 'red', label = "Mars' Orbit")

    ax.plot(0,0,0,marker='o',color='orange', label= 'Sun')
    ax.plot(R_EarthDP[0],R_EarthDP[1],R_EarthDP[2], marker='o', color='blue')
    ax.plot(R_MarsAR[0],R_MarsAR[1],R_MarsAR[2], marker='o', color='red')

    t, y = transfer_orbit.propagate_orbit()
    ax.plot(y[:, 0], y[:, 1], y[:, 2], color= 'magenta', label= 'Transfer Orbit')

    # Plot x,y,z axis vectors
    l = abs(y0[0] * 0.75)
    x,y,z = [[0,0,0], [0,0,0], [0,0,0]]
    u,v,w = [[l,0,0], [0,l,0], [0,0,l]]

    ax.quiver(x,y,z,u,v,w, color = 'lime', label= 'Inertial Frame Axes')

    ax.set_xlabel('X (10^8 km)')
    ax.set_ylabel('Y (10^8 km)')
    ax.set_zlabel('Z (10^8 km)')
    ax.set_zlim(-3e8, 3e8)
    plt.legend()
    plt.grid(color = 'lightgray', linestyle = '--', linewidth = 0.3)
    plt.show()


if __name__ == '__main__':
    main()
