# -*- coding: utf-8 -*-
"""
Created on Mon Nov 19 16:30:31 2018

@author: Dawn K Merriman
"""


import numpy as np
import mpmath as mp 

week_sec = 604800
WE = 7.2921151467e-5  # WGS 84 value of earth's rotation rate (rad/s)
frel = -4.442807633e-10 #Clock relativity parameter, (s/m^1/2)
mu = 3.986005e14 #Universal gravitational parameter (m^3/sec^2)

def GPSEph2xyz(df_epoch,gps_eph_df,tk_in):

    tk = (df_epoch.GPS_Week[0]-gps_eph_df.GPS_Week[0])*week_sec + (tk_in-df_epoch.Toe)

    df_epoch['A'] = df_epoch.Asqrt**2
    df_epoch['n0'] = np.sqrt(mu/(df_epoch.A**3))
    df_epoch['n'] = df_epoch.n0+df_epoch.delta_n
    df_epoch['h'] = np.sqrt(df_epoch.A*(1-df_epoch.e**2))*mu
    df_epoch['Mk'] = df_epoch.M0+df_epoch.n*tk  


    ### --- Kepler's Equation --- ### 
    df_epoch['err'] = 1    
    Ek = df_epoch.Mk

    count = 0
    max_count = 20
    while any(abs(df_epoch.err) > 1e-10) and count < max_count:
        df_epoch.err = Ek - df_epoch.Mk - df_epoch.e*\
            np.sin(Ek.astype(np.longdouble))
        Ek = Ek - df_epoch.err 
        count += 1
        #print count
        if count == max_count:
            print "Failed convergence on Kepler's equation"


    dt = (df_epoch.GPS_Week[0]-gps_eph_df.GPS_Week[0])*week_sec + (tk_in-df_epoch.Toc)

    #Switch to mpf format for better precision 
    sin_Ek = []
    cos_Ek = []
    for i in Ek:
        sin_Ek.append(mp.sin(i))
        cos_Ek.append(mp.cos(i))

    dtsvS = df_epoch.af0+df_epoch.af1*dt+df_epoch.af2*(dt**2)+frel\
        *df_epoch.e*df_epoch.Asqrt*sin_Ek-df_epoch.TGD
    
    #mp cannot deal with pandas types so need to convert all... BUH 
    vk = []
    for e,cos_ek,sin_ek,omega in zip(df_epoch.e,cos_Ek,sin_Ek,df_epoch.omega):
        y = mp.sqrt(1-mp.power(e,2))*sin_ek/(1-e*cos_ek)
        x = (cos_ek-e)/(1-e*cos_ek)
        vk.append(mp.atan2(y,x))
    
    Phik =  vk + df_epoch.omega

    duk = []
    drk = []
    dik = []
    for phik,cus,cuc,crc,crs,cis,cic in  zip(Phik,df_epoch.Cus,df_epoch.Cuc,df_epoch.Crc,df_epoch.crs,df_epoch.Cis,df_epoch.Cic):
        sin_2Phik = mp.sin(2*phik)
        cos_2Phik = mp.cos(2*phik)
        duk.append(cus*sin_2Phik+cuc*cos_2Phik) #Argument of latitude correction
        drk.append(crc*cos_2Phik+crs*sin_2Phik) #Radius Correction
        dik.append(cic*cos_2Phik+cis*sin_2Phik) #Correction to Inclination

    uk = Phik + duk #Corrected argument of latitude

    rk = []
    ik = []
    sin_uk = []
    cos_uk = []
    for A,e,VK,DRK,i,idot,TK,DIK,UK in zip(df_epoch.A,df_epoch.e,vk,drk,df_epoch.i0,df_epoch.IDOT,tk,dik,uk):
        rk.append(A*(1-mp.power(e,2))/(1+e*mp.cos(VK))+DRK) #Corrected radius
        ik.append(i+idot*TK+DIK) #Corrected inclination
        sin_uk.append(mp.sin(UK))
        cos_uk.append(mp.cos(UK))

    xkp = []
    ykp = []
    for RK,COS_UK,SIN_UK in zip(rk,cos_uk,sin_uk):
        xkp.append(RK*COS_UK) #Position in orbital plane
        ykp.append(RK*SIN_UK) #Position in orbital plane

    Wk = []
    for OMEGA,OMEGA_DOT,TK,TOE in zip(df_epoch.OMEGA,df_epoch.OMEGA_DOT,tk,df_epoch.Toe):
        Wk.append(OMEGA+(OMEGA_DOT - WE)*TK - WE*TOE) #Wk - corrected longitude of ascending node (WE is a constant)
    
    cartX = []
    cartY = []
    cartZ = []
    for wk,IK,XKP,YKP in zip(Wk,ik,xkp,ykp):
        cartX.append(XKP*mp.cos(wk)-YKP*mp.cos(IK)*mp.sin(wk))
        cartY.append(XKP*mp.sin(wk)+YKP*mp.cos(IK)*mp.cos(wk))
        cartZ.append(YKP*mp.sin(IK))
    
    return cartX,cartY,cartZ,dtsvS