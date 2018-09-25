# -*- coding: utf-8 -*-
"""
Created on Thu Sep 06 11:08:30 2018

@author: Dawn K Merriman

"""
from Read_GNSS_log_v1 import *
from Reading_RINEX_Nav_v3 import *
from math import *
import pandas as pd
import numpy as np
import mpmath as mp
import matplotlib.pyplot as plt
from datetime import datetime 
from matplotlib.backends.backend_pdf import PdfPages

startTime = datetime.now()

#As of 1 Jan 2017 (http://tycho.usno.navy.mil/leapsec.html) 
#GPS-UTC = 18 seconds 
lsec = 18

data_file = 'C:\Users\Cameron\Dropbox\PostDoc_Work\Thule_Files\Thule\gnss_log_2018_06_01_09_43_31_dkm.txt'
#gps_rinex_file = 'C:\Users\dawnk\Dropbox\PostDoc_Work\Thule_Files\Ephemeris_data\hour1520.18n'
#glonass_rinex_file ='~\Dropbox\PostDoc_Work\Thule_Files\Ephemeris_data\hour1520.18g'
#other_rinex_file = 'C:\Users\Cameron\Dropbox\PostDoc_Work\Thule_Files\'

mp.dps = 64

## -- Constants WGA 84 fundamental parameters revised in 1997 -- ##
a = 6378137.0 # semi-major axis of ellipsoid in m
f = 298.257223563 # 1/f the reciprocal flattening 
omga_e = 7292115e-11 #Earth's angular velocity in rad/sec
grav = 3986004.418e-8 #Earth's gravational constant m^3/s^2
c = 2.99792458e8 #speed of light m/s

##Time Keeper Constants##
min_sec = 60
hour_sec = 3600
day_sec = 86400
week_sec = 604800
month_days = [31,28,31,30,31,30,31,31,30,31,30,31] #not leap year months
lsec = np.float(18.0/60)

light_speed = 299792458 #m/s


df = readGNSSLogger(data_file) #Read in file
#gps_eph_df = get_eph(gps_rinex_file) #Read in ephmersis data for GPS
#glonass_eph_df = get_eph(glonass_rinex_file) #Glonass Ephmersis data

#### ---------------------- FILTERING DATA ----------------- ########
df = df[df.ConstellationType==1] #SELECTING ONLY GPS DATA
#df = df[df.Svid == 31] #SELECTING ONLY ONE SATELLITE AT THE MOMENT#

### REMOVING ALL DATA USING BIT NOTATION BASED ON THE STATE OF THE GNSS ENGINE ###
df['my_index'] = df.index
bad_data_rows = [] 
for index,row in df.iterrows():
    #print index ### NOTE IF USING PANDAS IT WILL ORDER BASED ON SMALLEST NUMBER NOT ON INDEX UNLESS YOU SPECIFIY 
    if int(row['State']) & 2**3 != 8 and int(row['State']) & 2**0 != 1:
        bad_data_rows.append(row['my_index'])
    #elif abs(row['PseudorangeRateUncertaintyMetersPerSecond']*2) > row['PseudorangeRateMetersPerSecond']:
     #   bad_data_rows.append(row['my_index'])

df = df.drop(bad_data_rows)

datafiltertime = datetime.now()-startTime

print 'Done Filtering Data'
print datafiltertime
## HOLY SHIT IT TAKES A LONG TIME TO FILTER THIS DATA ## 

##GPS to UTC time 
#Convert GPS time to [yyyy,mm,dd,hh,mm,ss] with no leap secs
df['allRxMillis'] = (df.TimeNanos - df.FullBiasNanos)* pow(10,-6) #True GPS time ie time since 1980/1/1 until now in ms
df['days'] = (df.allRxMillis*pow(10,-3)/day_sec) + 5 #days since 1980/1/1

leap_count = []
year = 1980
c = 1

#days = df.days[0] #taking time of 1st data point
days = df.days.iloc[0] #taking the time of last data point

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
tmin = np.floor(((days-day)*24-thour)*60-lsec)
tsec = round(((((days-day)*24-thour)*60-lsec)-tmin)*60,1)
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
df['PrM_sigma'] = df.ReceivedSvTimeUncertaintyNanos*light_speed

df['ElapsedTimehr'] = (df.ElapsedRealtimeMillis-min(df.ElapsedRealtimeMillis))*10**-3/60/60

### REMOVING PRM = 0 DATA ### 
for index,row in df.iterrows():
    if int(row['PrM']) == 0:
        df.drop(index, inplace=True)

datafiltertime = datetime.now()-startTime

print 'Done Filtering Data v2'
print datafiltertime

## -- MORE DATA FILTERING -- ## 
bad_data_rows = [] 
for index,row in df.iterrows():
    #print index ### NOTE IF USING PANDAS IT WILL ORDER BASED ON SMALLEST NUMBER NOT ON INDEX UNLESS YOU SPECIFIY 
    if int(row['PrM']) > 1e10: ## again more data that looks terrible ## 
        bad_data_rows.append(index)

df = df.drop(bad_data_rows)

uniquesvids = np.unique(df.Svid)
print uniquesvids 

sat_names=["GPS","SBAS","GLONASS","QZSS","BEIDOU","GALILEO","UNKNOWN","UNKNOWN","UNKNOWN","UNKNOWN"]

with PdfPages('gnss_log_2018_06_04_13_37_47_pseudoranges.pdf') as pdf:
    for svid in uniquesvids:
        subdf = df[df.Svid==svid]
        
        plt.figure()
        fig, ax = plt.subplots()
        plt.title("{0} #{1}".format(sat_names[int(subdf.ConstellationType.iloc[0])-1],svid))
        plt.plot(subdf.ElapsedTimehr,subdf.PrM)
        plt.xlabel('Elapsed Real Time (hour)')
        plt.ylabel('Pseudorange (m)')
        pdf.savefig()
        plt.close
        
        #plt.legend(uniquesvids)

plt.figure()
plt.xlim(0,1.6)
for svid in uniquesvids:
    subdf = df[df.Svid==svid]
    
    plt.plot(subdf.ElapsedTimehr,subdf.PrM)
    plt.xlabel('Elapsed Real Time (hour)')
    plt.ylabel('Pseudorange (m)')
    #plt.legend(uniquesvids)


pdf.savefig()
plt.close
fig.savefig('gnss_log_2018_06_04_13_37_47_all_pseduoranges.pdf', \
bbox_inches='tight')
            
        
        
    
