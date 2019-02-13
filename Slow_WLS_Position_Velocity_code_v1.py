# -*- coding: utf-8 -*-
"""
Created on Thu Sep 06 11:08:30 2018
 
@author: Dawn K Merriman
 
"""
from Read_GNSS_log_v1 import *
from Reading_RINEX_Nav_v3 import * 
from Find_Closest_Epoch import *
from GPS_Eph_to_XYZ_v3 import *
from FlightTimeCorrection import *
from RotEcef2Ned import *
from Xyz2Lla import *
from numpy import linalg as LA
from math import *
import pandas as pd
import numpy as np
import mpmath as mp
from coordinates import * #this is taken from AENeAS code

import matplotlib.pyplot as plt
from datetime import datetime 
from matplotlib.backends.backend_pdf import PdfPages
 
startTime = datetime.now()
 
print 'Note: Leap second is currently set to 18 secs'
 
#data_file = '/home/dawn/googledrive/Work/GNSS_AFIT_Project/Raw_Data_Files/Thule/gnss_log_2018_06_01_09_43_31_dkm.txt'
#gps_rinex_file = '/home/dawn/googledrive/Work/GNSS_AFIT_Project/Raw_Data_Files/Thule/hour1520.18n'
 
data_file = 'D:\Dawn\Google_Drive\Work\GNSS_AFIT_Project\Raw_Data_Files\Thule\gnss_log_2018_06_01_09_44_05_dkm_trunc.txt'
gps_rinex_file = 'C:\Users\Cameron\Dropbox\PostDoc_Work\Thule_Files\Ephemeris_data\hour1520.18n'

#glonass_rinex_file ='~\Dropbox\PostDoc_Work\Thule_Files\Ephemeris_data\hour1520.18g'
#other_rinex_file = 'C:\Users\Cameron\Dropbox\PostDoc_Work\Thule_Files\BRDC00WRD_R_20181520000_01D_CN.rnx'
 
mp.dps = 64
 
## -- Constants WGA 84 fundamental parameters revised in 1997 -- ##
f = 298.257223563 # 1/f the reciprocal flattening 
omga_e = 7292115e-11 #Earth's angular velocity in rad/sec
grav = 3986004.418e-8 #Earth's gravational constant m^3/s^2
 
mu = 3.986005e14 #Universal gravitational parameter (m^3/sec^2)
eartheccen2 = 6.69437999014e-3 #Earth eccentricity squared m^2
earthmeanradius = 6371009 #Mean R of ellipsoid(m) IU Gedosey& Geophysics
earthsemimajor = 6378137 #Earth semi-major axis (m)
ephvalidsec = 7200 #+- 2 hours ephemeris validity
frel = -4.442807633e-10 #Clock relativity parameter, (s/m^1/2)
gpsepochjd = 2444244.5 #GPS Epoch in Julian Days
horizendeg = 5 #angle above horizon at which GPS models break down
meanflightsecs = 75e-3 #mean time of flight btwn closest GPS sat (~66 ms) & furthest (~84 ms):
WE = 7.2921151467e-5  # WGS 84 value of earth's rotation rate (rad/s)
 
##Time Keeper Constants##
min_sec = 60
hour_sec = 3600
day_sec = 86400
week_sec = 604800
month_days = [31,28,31,30,31,30,31,31,30,31,30,31] #not leap year months
 
 
#As of 1 Jan 2017 (http://tycho.usno.navy.mil/leapsec.html) 
#GPS-UTC = 18 seconds 
lsec = 18
lsec_min= np.float(18/60)
 
light_speed = 299792458 #m/s 

#Thresholds
maxdelposfornavm = 20
maxprruncmps = 10
maxtowuncns = 500
mingpseph = 24
 
 
df = readGNSSLogger(data_file) #Read in file
gps_eph_df = get_eph(gps_rinex_file) #Read in ephmersis data for GPS
#glonass_eph_df = get_eph(glonass_rinex_file) #Glonass Ephmersis data
 
print np.unique(df.ConstellationType)
 
#### ---------------------- FILTERING DATA ----------------- ########
#df = df[df.ConstellationType==1] #SELECTING ONLY GPS DATA
#df = df[df.Svid == 31] #SELECTING ONLY ONE SATELLITE AT THE MOMENT#
 
### REMOVING ALL DATA USING BIT NOTATION BASED ON THE STATE OF THE GNSS ENGINE ###
df['my_index'] = df.index
#bad_data_rows = [] 
#for index,row in df.iterrows():
    #print index ### NOTE IF USING PANDAS IT WILL ORDER BASED ON SMALLEST NUMBER NOT ON INDEX UNLESS YOU SPECIFIY 
#    if int(row['State']) & 2**3 != 8 and int(row['State']) & 2**0 != 1:
#        bad_data_rows.append(row['my_index'])
    #elif abs(row['PseudorangeRateUncertaintyMetersPerSecond']*2) > row['PseudorangeRateMetersPerSecond']:
        #bad_data_rows.append(row['my_index'])
 
#df = df.drop(bad_data_rows)
print np.unique(df.ConstellationType)
 
datafiltertime = datetime.now()-startTime
 
print 'Done Filtering Data'
print datafiltertime
 
##GPS to UTC time 
#Convert GPS time to [yyyy,mm,dd,hh,mm,ss] with no leap secs
df['allRxMillis'] = (df.TimeNanos - df.FullBiasNanos)* pow(10,-6) #True GPS time ie time since 1980/1/1 until now in ms
df['days'] = (df.allRxMillis*pow(10,-3)/day_sec) + 5 #days since 1980/1/1
 
leap_count = []
year = 1980
c = 1
 
#days = df.days[0] #taking time of 1st data point
days = df.days.iloc[0] #taking the time of first data point
 
while days >= 366:
    if c < 4:
        days = days - 365
        c += 1
        leap_count.append(0)
    elif c == 4:
        days = days - 366
        c = 1
        leap_count.append(1)
year = np.floor(1980 + np.shape(leap_count)[0])
 
if leap_count[-1] == 1: month_days[1] = 29 #change days in Feb is current year is a leap year
 
month = 0
while days >= month_days[month]:
    days = days - month_days[month]
    month += 1
month += 1
day = np.floor(days)
thour = np.floor((days-day)*24)
tmin = np.floor(((days-day)*24-thour)*60-lsec_min)
tsec = round(((((days-day)*24-thour)*60-lsec_min)-tmin)*60,1)
gps_week = np.floor(-df.FullBiasNanos*pow(10,-9)/week_sec)
 
print year,month,day,thour,tmin,tsec
 
##--- Processing ----### 
 
## -- Not sure how the State reflects if the TOW and FullBiasNanos is correct ?? What are the State numbers? -- ##
 
#check that df.FullBiasNanos is correct by checking if df.State[0] has bits 0 and 3 if not find first datapoint that does
for i in df.index: ### NOTE IF USING PANDAS IT WILL ORDER BASED ON SMALLEST NUMBER NOT ON INDEX UNLESS YOU SPECIFIY 
    if df.State[i] & 2**3:
        FullBiasNanos = df.FullBiasNanos[i]
        break
    else:
        pass
     
df['tRxNanos'] = df.TimeNanos - FullBiasNanos - (gps_week*mp.mpf(week_sec*pow(10,9)))
 
#check if any tRxNanos are negative#
if any(df.tRxNanos < 0) == True: print 'Negative time detected! ...HELP!'
 
#subtract the fractional offsets TimeOffsetNanos and BiasNanos
df['tRxSecs'] = (df.tRxNanos - df.TimeOffsetNanos - df.BiasNanos)*1e-9
df['tTxSecs'] = df.ReceivedSvTimeNanos*1e-9
 
##--- Check for week roll over ---##
df['prSeconds'] = df.tRxSecs - df.tTxSecs
 
#if df.prSeconds > week_sec: print 'SHIT THERE IS A WEEK ROLLOVER... PLS CORRECT ME'
 
## Calc Pseudoranges ## 
df['PrM'] = df.prSeconds*light_speed
df['PrM_sigma'] = df.ReceivedSvTimeUncertaintyNanos*light_speed*1e-9
 
df['ElapsedTimehr'] = (df.ElapsedRealtimeMillis-min(df.ElapsedRealtimeMillis))*10**-3/60/60
df['tS'] = df.tRxSecs - df.prSeconds
 
### REMOVING PRM = 0 DATA ### 
for index,row in df.iterrows():
    if int(row['PrM']) == 0:
        df.drop(index, inplace=True)
 
datafiltertime = datetime.now()-startTime
 
print 'Done Filtering Data v2'
print datafiltertime
print np.unique(df.ConstellationType)
 
## -- MORE DATA FILTERING -- ## 
bad_data_rows = [] 
for index,row in df.iterrows():
    #print index ### NOTE IF USING PANDAS IT WILL ORDER BASED ON SMALLEST NUMBER NOT ON INDEX UNLESS YOU SPECIFIY 
    if int(row['PrM']) > 1e14: ## again more data that looks terrible ## 
        bad_data_rows.append(index)
 
df = df.drop(bad_data_rows)
print 'Done Filtering Data v3'
 
datafiltertime = datetime.now()-startTime
print np.unique(df.ConstellationType)
print datafiltertime
 
uniqueconstellations = np.unique(df.ConstellationType)
 
sat_names=["GPS","SBAS","GLONASS","QZSS","BEIDOU","GALILEO","UNKNOWN","UNKNOWN","UNKNOWN","UNKNOWN"]
 
#with PdfPages('test_pseudoranges.pdf') as pdf:
#    for cons in uniqueconstellations:
#        print cons
#        print sat_names[int(cons)-1]
#        df2 = df[df.ConstellationType==cons] #SELECTING ONLY GPS DATA
#        uniquesvids = np.unique(df.Svid)
#        #print uniquesvids 
#        for svid in uniquesvids:
#            subdf = df2[df2.Svid==svid]
#         
#            plt.figure()
#            fig, ax = plt.subplots()
#            plt.title("{0} #{1}".format(sat_names[int(cons)-1],svid))
#            plt.plot(subdf.ElapsedTimehr,subdf.PrM)
#            plt.xlabel('Elapsed Real Time (hour)')
#            plt.ylabel('Pseudorange (m)')
#            pdf.savefig()
#            plt.close
#         
#            plt.legend(uniquesvids)
# 
#plt.figure()
#plt.xlim(0,1.6)
#for cons in uniqueconstellations:
#    df2 = df[df.ConstellationType==cons] #SELECTING ONLY GPS DATA
#    uniquesvids = np.unique(df.Svid)
#    for svid in uniquesvids:
#        subdf = df[df.Svid==svid]
#     
#        plt.plot(subdf.ElapsedTimehr,subdf.PrM, label = "{0} #{1}".format(sat_names[int(cons)-1],svid))
#        plt.xlabel('Elapsed Real Time (hour)')
#        plt.ylabel('Pseudorange (m)')
#        plt.legend()
# 
# 
#pdf.savefig()
#plt.close
#fig.savefig('test_pseduoranges.pdf', \
#bbox_inches='tight')
# 
#### -------------------------- Eph data analysis --------------------------###
if any(gps_eph_df.year) > 80:
    gps_eph_df.year = 1900+gps_eph_df.year
else:
    gps_eph_df.year = 2000+gps_eph_df.year
     
### ---------- UTC to GPS time -------------###   
gps_eph_df['JulianDay'] = np.floor(365.25*gps_eph_df.year.astype(np.float64)) + \
    np.floor(30.6001*(gps_eph_df.month.astype(np.float64)+1)) - 15 + 1720996.5 + \
    gps_eph_df.day.astype(np.float64) + gps_eph_df.hour.astype(np.float64)/24
gps_eph_df['DaysSinceEpoch'] = np.floor(gps_eph_df.JulianDay.astype(np.float64)\
          - gpsepochjd)
gps_eph_df['gpsweek'] = (gps_eph_df.DaysSinceEpoch/7)
gps_eph_df['dayofweek'] = np.remainder(gps_eph_df.DaysSinceEpoch,7)
gps_eph_df['Toc'] = gps_eph_df.dayofweek*day_sec + gps_eph_df.hour*hour_sec + \
    gps_eph_df.minutes*min_sec + gps_eph_df.sec+lsec
     
# get time from measurements 
FctSeconds = df.allRxMillis*10**-3
FctSeconds = np.unique(FctSeconds)
 
##-- Check what ephemeris data is needed --## 
#for sat in np.unique(df.ConstellationType):
#    if sat == 1:
         
## --- Get SVIDs for XX CONSTELLATION --- ## NEED TO CYCLE CONSTELLATIONS HERE
time_start = datetime.now()
df_gps = df[df.ConstellationType==1] ##only grabbing GPS at the moment
gps_svids = np.unique(df_gps.Svid)
if np.size(gps_svids) < 4:
    print 'Shit, there are not enough GPS satelites to calculate a poisiton!'
    exit()
     
### NEED TO CYCLE EACH UNIQUE FCTSECOND STARTING HERE ###

#Variabls we are saving from each FctSecond cycled 
allLlaDegDegM = []
allBcMeters = []
allVelMps = []
allBcDotMps = []
hdop = []
sigmaLLaM = []
sigmaVelMps = []

start_slove_position = datetime.now()

datafiltertime = datetime.now()-startTime
print datafiltertime

points = np.shape(FctSeconds)[0]
info_count = 0

for time in FctSeconds:
    
    info_count += 1
    
    #print "You are on point {0} out of {1}".format(info_count,points)
    
    ## Get closest epoch bassed on FctSecond measurement ##
    min_index = close_eph(gps_svids,gps_eph_df,time)
 
    #creating new dataframes for this time measurement for measured and ephemeris data
    df_fct_gps = df_gps[df_gps.allRxMillis*10**-3 == time]
    first_epoch = gps_eph_df.loc[min_index] 
    #find difference between sets of Svids
    dif = list(set(first_epoch.svid) - set(df_fct_gps.Svid)) 
    #take out GPS Svids that do not have measurements as well 
    for i in dif:
        test_df = first_epoch[first_epoch.svid == i]
        first_epoch = first_epoch.drop(test_df.index)

    if np.size(first_epoch.index) < 4:
        print "Shit, first epoch doesn't have enough satsllites to calc position"
        exit()
    ## CHANGE THIS TO GO TO NEXT FCTSECOND ## 

    #supress pandas SettingWithCopyWarning -- no idea why it pops up for the simple
    #subtraction below... 
    pd.options.mode.chained_assignment = None  # default='warn'
 
    df_fct_gps['ttxSeconds'] = df_fct_gps.tRxSecs - df_fct_gps.PrM/light_speed #fixed 1/30/19
    df_fct_gps = df_fct_gps.reset_index(drop=True) 
    first_epoch = first_epoch.reset_index(drop=True) 

    ### GpsEph2Dtsv.m 
    ### ---- Calculate dependent variables ---- ###
    gpsWeek = np.floor((df_fct_gps.allRxMillis*10**-3)/week_sec)
    tk = df_fct_gps.ttxSeconds - first_epoch.Toe 
 
    #time since time of applicability 
    ##Check if there is roll over ??##
    #302400 = week_sec/2
    cor_tk = []
    if any(abs(tk) > 302400) == True:
        for time in tk:
            if time > 302400:
                cor_tk.append(time - week_sec)
            if time < 302400:
                cor_tk.append(time + week_sec)
    else:
        cor_tk = tk

    first_epoch['A'] = first_epoch.Asqrt**2
    first_epoch['n0'] = np.sqrt(mu/(first_epoch.A**3))
    first_epoch['n'] = first_epoch.n0+first_epoch.delta_n
    first_epoch['h'] = np.sqrt(first_epoch.A*(1-first_epoch.e**2))*mu
    first_epoch['Mk'] = first_epoch.M0+first_epoch.n*cor_tk  
    #first_epoch['Ek'] = first_epoch.Mk
 
    ### --- Kepler's Equation --- ### 
    first_epoch['err'] = 1
    Ek = first_epoch.Mk
    count = 0
    max_count = 20
    while any(abs(first_epoch.err) > 1e-8) and count < max_count:
        first_epoch.err = Ek - first_epoch.Mk - first_epoch.e*\
            np.sin(Ek.astype(np.float64))
        Ek = Ek - first_epoch.err 
        count += 1
        #print count
        if count == max_count:
            print "Failed convergence on Kepler's equation"
 
    dt = df_fct_gps.ttxSeconds - first_epoch.Toc
    cor_dt = []
    if any(abs(dt) > 302400) == True:
        for time in dt:
            if time > 302400:
                cor_dt.append(time - week_sec)
            if time < 302400:
                cor_dt.append(time + week_sec)
    else:
        cor_dt = dt

    #clock bias
    dtsvS=first_epoch.af0+first_epoch.af1*cor_dt+first_epoch.af2*(cor_dt**2)+frel\
        *first_epoch.e*first_epoch.Asqrt*np.sin(Ek.astype(np.float64)) - first_epoch.TGD
    
    df_fct_gps['ttx'] = df_fct_gps.ttxSeconds - dtsvS
    df_fct_gps['week'] = np.floor(-df_fct_gps.FullBiasNanos*pow(10,-9)/week_sec)

    ##---- GpsEph2Pvt.m  ------###
    #compute velocity from delta position and dtsvS at (t+0.5) - (t-0.5)
    #This is better than differentiating beacuse both the orbital and relativity
    #terms have nonlinearities that are not easily differentiable 
    x,y,z,dtsvS = GPSEph2xyz(first_epoch,gps_eph_df,df_fct_gps.ttx)
    t1 = list(df_fct_gps.ttx + 0.5) 
    xp,yp,zp,dtsvSp = GPSEph2xyz(first_epoch,gps_eph_df,t1)
    t2 = list(df_fct_gps.ttx - 0.5)
    xm,ym,zm,dtsvSm = GPSEph2xyz(first_epoch,gps_eph_df,t2)

    vMpsx = []
    vMpsy = []
    vMpsz = []
    dtsvSDot = []
    for XP,YP,ZP,DTSVSP,XM,YM,ZM,DTSVSM in zip(xp,yp,zp,dtsvSp,xm,ym,zm,dtsvSm):
        vMpsx.append(XP-XM)
        vMpsy.append(YP-YM)
        vMpsz.append(ZP-ZM)
        dtsvSDot.append(((DTSVSP-DTSVSM)*light_speed)) #NOTE YOU ARE MULTIPLYING BY speed of light here to avoid memory err later on in code 

    ##---- back to WlsPvt.m ------###
    #Compute weights 
    Wpr = np.diag(1/df_fct_gps.PrM_sigma)
    Wrr = np.diag(1/df_fct_gps.PseudorangeRateUncertaintyMetersPerSecond)

    s = np.size(df_fct_gps.PseudorangeRateMetersPerSecond)

    bc = 0
    xyz0 = np.matrix('0;0;0')
    xo = np.zeros((8,1))
    xHat =  np.zeros((4,1))
    dx = xHat+float("inf")

    count = 0
    max_count = 20
    while LA.norm(dx) > maxdelposfornavm:
        count += 1

        svXyzTrx = []
        dtflight_Test = []
        for i in range(0,s):
            dtflight_Test.append((df_fct_gps.PrM[i] - bc)/light_speed+dtsvS[i])
            dtflight = (df_fct_gps.PrM[i] - bc)/light_speed+dtsvS[i]
            svXyzTrx.append(FlightTimeCorrection(x[i],y[i],z[i],dtflight))

        svXyzTrx = np.squeeze(svXyzTrx, axis=0) #getting rid of list dimension needed for loop 
        v = xyz0*np.ones(s) - np.matrix.transpose(svXyzTrx)
        v = np.matrix.transpose(v)
        range_v = np.sqrt(np.matrix.sum(np.square(v),axis=1))
        v = np.divide(v,range_v) #line of sight unit vectors from sv to xo

        prHat = np.squeeze(range_v) + bc - light_speed*dtsvS.values
        
        zPr = df_fct_gps.PrM.values - prHat
        H = np.hstack((v,np.ones((s,1)))) #H matirx [unit vector,1]

        #dx is off from what MatLab says in WlsPvt line 110 
        dx = np.linalg.pinv((Wpr*H).astype(np.float64))*Wpr*np.transpose(zPr) #z = Hx, premultiply by W: Wz = WHx, and solve for x

        #update values
        xHat = xHat + dx
        xyz0 = xyz0 + dx[0:3]
        bc = bc + np.squeeze(np.array(dx[3]))

        #Calculate a-posteriori range residual 
        zPr = np.transpose(zPr)-H*dx

    #Compute velocities --------- 
    rrMps = []
    for i in range(0,s):
        rrMps.append(-1*np.matrix([vMpsx[i],vMpsy[i],vMpsz[i]])*np.transpose(v[i]))

    prrHat = np.squeeze(rrMps) + xo[7] - dtsvSDot
    zPrr = df_fct_gps.PseudorangeRateMetersPerSecond - prrHat
    vHat = np.linalg.pinv((Wrr*H).astype(np.float64))*Wrr*np.transpose(np.asmatrix(zPrr.values))
    xHat = np.concatenate((xHat,vHat),axis=0)
    z = np.concatenate((zPr,np.transpose(np.asmatrix(zPrr.values))),axis=0)

    #Back to GpsWlsPvt.m line 85 
    xo = (xo + xHat).astype(np.float64)

    #extract position states 
    llaDegDegM = np.squeeze(Xyz2Lla(xo[0],xo[1],xo[2]))
    allLlaDegDegM.append(llaDegDegM)
    allBcMeters.append(xo[3])

    #extract velocity states 
    RE2N = RotEcef2Ned(llaDegDegM[0],llaDegDegM[1])
    vNed = RE2N*xo[4:7]
    allVelMps.append(vNed)
    allBcDotMps.append(xo[7])

    #compute HDOP
    H = np.hstack((H[:,0:3]*np.transpose(RE2N),np.ones((s,1))))
    P = np.linalg.inv((np.transpose(H)*H).astype(np.float64))
    hdop.append(np.sqrt(P[0,0]+P[1,1]))

    #Compute variance of llaDegDegM
    P = np.linalg.inv((np.transpose(H)*(np.transpose(Wpr)*Wpr)*H).astype(np.float64))
    sigmaLLaM.append(np.sqrt(np.diag(P[0:3,0:3])))

    #Computer variance of velocity 
    P = np.linalg.inv((np.transpose(H)*(np.transpose(Wrr)*Wrr)*H).astype(np.float64))
    sigmaVelMps.append(np.sqrt(np.diag(P[0:3,0:3])))

total_time = datetime.now()-startTime
position_time = datetime.now()-start_slove_position
datafilter_time = datafiltertime

print "Total time = {0} \n Time to filter data = {1} \n Time to calc position = {2}".format(total_time,datafilter_time,position_time)
