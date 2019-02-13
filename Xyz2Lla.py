# -*- coding: utf-8 -*-
"""
Created on Tue Feb 12 12:37:55 2019

@author: Dawn K Merriman
"""

import numpy as np

def Xyz2Lla(x, y, z):
    
    if x == 0 or y == 0:
        print "x or y is zero, lla is undefined" 
        llaDegDegM = np.array([np.nan,np.nan,np.nan])
    
    e2 =  6.69437999014e-3 #Earth eccentricity squared m^2
    a = 6378137 #Earth semi-major axis (m)
    a2 = pow(6378137,2)
    b2 = a2*(1-e2)
    b = np.sqrt(b2)
    ep2 = (a2 - b2)/b2
    p = np.sqrt(x**2+y**2)
    
    s1 = z*a
    s2 = p*b
    h = np.sqrt(s1**2+s2**2)
    sin_theta = s1/h
    cos_theta = s2/h
    
    s1 = z+ep2*b*(sin_theta**3)
    s2 = p-a*e2*(cos_theta**3)
    h = np.sqrt(s1**2 + s2**2)
    tan_lat = s1/s2
    sin_lat = s1/h
    cos_lat = s2/h
    latDeg = np.arctan(tan_lat)
    latDeg = latDeg*(180/np.pi)
    

    N = a2*(1/np.sqrt(a2*(cos_lat**2) + b2*(sin_lat**2)))
    altM = p/cos_lat - N
    
    lonDeg = np.arctan2(y,x)
    lonDeg = np.remainder(lonDeg,2*np.pi)*(180/np.pi)
    if lonDeg > 180:
        lonDeg = lonDeg - 360
    
    llaDegDegM = np.array([latDeg, lonDeg, altM])
    
    
    return llaDegDegM