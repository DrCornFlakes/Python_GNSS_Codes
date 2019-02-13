# -*- coding: utf-8 -*-
"""
Created on Tue Feb 12 12:07:11 2019

@author: Dawn K Merriman 
"""

import numpy as np

def RotEcef2Ned(latDeg, lonDeg):
    
    D2R = np.pi/180;
    latRad = D2R*latDeg
    lonRad = D2R*lonDeg
    
    clat = np.cos(latRad)
    slat = np.sin(latRad)
    clon = np.cos(lonRad)
    slon = np.sin(lonRad)
    
    Re2n = np.array([[-slat*clon, -slat*slon, clat], [-slon, clon, 0], [-clat*clon, -clat*slon, -slat]])
    
    return Re2n