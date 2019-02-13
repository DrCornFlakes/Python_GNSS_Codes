# -*- coding: utf-8 -*-
"""
Created on Fri Feb  1 14:41:53 2019

@author: Dawn K Merriman
"""

##--- Flight time corrections ---##
import numpy as np
import mpmath as mp


WE = 7.2921151467e-5  # WGS 84 value of earth's rotation rate (rad/s)

def FlightTimeCorrection(x,y,z,dTflightSeconds):
    
    xE = np.matrix([x,y,z]).astype(np.float64)
    
    theta = WE*dTflightSeconds

    #Rotation matrix for ECEF to ECI 
    #opposite of what is in GPS 200-E 
    R3 = np.matrix([[mp.cos(theta), mp.sin(theta), 0], \
                [-1*mp.sin(theta), mp.cos(theta), 0], \
                [0, 0, 1]]).astype(np.float64)

    xERot = np.matmul(xE,R3)
    
    return xERot