#---------- Gregorian Date Calculator ----------#

# Calculation to find gregorian date given a juliandate. Using formula provided by University of Texas at Austin:
#   https://quasar.as.utexas.edu/BillInfo/JulianDatesG.html

import math as m
import numpy as np

def gregorian(JD):

    Q = JD + 0.5
    Z = Q//1
    W = int((Z - 1867216.25)/36524.25)
    X = int(W/4)
    A = Z + 1 + W - X
    B = A + 1524
    C = int((B - 122.1)/365.25)
    D = int(365.25*C)
    E = int((B-D)/30.6001)
    F = int(30.6001*E)
    Day = int(B-D-F+(Q-Z))

    # Must get number less than or equal to 12 for month
    if E-1 <= 12:
        Month = E-1
    else:
        Month = E-13

    # Year = C-4715 (if Month is January or February) or C-4716 (otherwise)
    if Month <= 2:
        Year = C-4715
    else:
        Year = C-4716

    

    return Month,Day,Year
