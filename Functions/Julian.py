#---------- Julian Date Calculator ----------#


import math as m
import numpy as np

def julianpre(year,month,day,hour,min):
    JD = 367*year-m.floor((7*(year+m.floor((month+9)/12)))/4)+m.floor(275*month/9)+day+1721013.5+(((5/60+min)/60)+hour)/24
    format(JD,'.17f')
    return JD

def julianinp():
    year = float(input("Enter the 4-digit year you are looking to convert: "))
    month = float(input("Enter the month of the year you are looking to convert (1-12): "))
    day = float(input("Enter the day of conversion: "))
    hour = float(input("Enter the hour of the day of the conversion (in military time 1-23): "))
    min = float(input("Enter how many minutes into the hour (1-59): "))

    JD = 367*year-m.floor((7*(year+m.floor((month+9)/12)))/4)+m.floor(275*month/9)+day+1721013.5+(((5/60+min)/60)+hour)/24
    format(JD,'.17f')
    return JD
