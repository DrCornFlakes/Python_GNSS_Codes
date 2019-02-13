# -*- coding: utf-8 -*-
"""
Created on Mon Oct 15 19:03:46 2018

@author: Dawn K Merriman
"""
import pandas as pd

## --- Finding the cloest epoch in eph data to measured data --- ##
def close_eph(svids,eph_data_frame,FctSeconds):
    week_sec = 604800

    min_index = []
    for sat in svids:
        eph_pull = eph_data_frame[eph_data_frame.svid==sat] 
        fctToe = eph_pull.GPS_Week*week_sec + eph_pull.Toe
        subtract_matrix = abs(fctToe - FctSeconds)
        subtract_matrix = subtract_matrix.astype('float')
        check = subtract_matrix.idxmin()
        if min(subtract_matrix) < (eph_data_frame.Fit_interval.loc[check]/2)*3600 \
            or (eph_data_frame.Fit_interval.loc[check] == 0 and min(subtract_matrix) < 2*3600):
            min_index.append(subtract_matrix.idxmin()) #pull index of min value 
        else:
            print 'No valid ephemeris found for svid %s' % sat
    return min_index

#Note that Rinex states the fit intervale == 0 it is unknown so we set it to 
#the typical 4 hour interval (see line 22)

#### --- This section of code copies what the MatLab code does may not be best
#### --- for our code 
## creating a matrix of each Toe time and measured time ##
#subtract_matrix = abs(fctToe - FctSeconds) #using only first measured time 
#need to have a list to find min easily#
#Note to self:
#This could be made faster but only using the first measured data point instead of all 5,600
#is that better? --- for us maybe not becuase we will later be doing 'real-time'?
#flat_list = [item for sublist in subtract_matrix.values.tolist() for 
#             item in sublist]
##getting the index number so we can pull from the eph file## 
#for index in list(subtract_matrix.index):
#    if any(subtract_matrix.loc[index] == min_value):
#        min_index = index

