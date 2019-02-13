# -*- coding: utf-8 -*-
"""
Created on Tue Sep 04 14:25:18 2018

@author: Dawn K Merriman
"""

import mpmath as mp
import pandas as pd

##### ----- Reading RINEX Navigation Broadcast ephemeris data ------ #####

#D:\Dawn\PostDoc_Work\Thule_Files\Ephemeris_Data\hour1520.18n

def get_eph(file_path):
    
    eph_h = {}
    eph=[]
    eph1=[]
    eph2=[]
    eph3=[]
    eph4=[]
    eph5=[]
    eph6=[]
    eph7=[]

    with open(file_path,'r') as f:
        cnt = 1
        for i in f:
            i = i.replace('D','E')
            if 'ION ALPHA' in i:
                eph_h['ion_alpha'] = mp.mpf(i[3:14]),mp.mpf(i[15:26]),mp.mpf(i[27:38]),mp.mpf(i[39:50])
            elif 'ION BETA' in i:
                eph_h['ion_beta'] = mp.mpf(i[3:14]),mp.mpf(i[15:26]),mp.mpf(i[27:38]),mp.mpf(i[39:50])
            elif 'DELTA-UTC' in i:
                eph_h['delta-utc'] = mp.mpf(i[3:22]),mp.mpf(i[23:40]), mp.mpf(i[44:50]), mp.mpf(i[55:59])
            elif 'LEAP SECONDS' in i:
                eph_h['leap_sec'] = int(i[3:30])
        ##--- Done with Header ---##
            elif any(char.isdigit() for char in i[0:2]):
                eph.append((mp.mpf(i[0:2]),mp.mpf(i[2:5]),mp.mpf(i[6:8]),mp.mpf(i[9:11]),mp.mpf(i[12:14]),mp.mpf(i[15:18]),mp.mpf(i[19:22]),mp.mpf(i[22:41]),mp.mpf(i[41:60]),mp.mpf(i[60:80])))  
            elif i[0:2]=='  ' and not 'A' in i:
                if cnt == 1:
                    eph1.append((mp.mpf(i[3:22]),mp.mpf(i[22:41]),mp.mpf(i[41:60]),mp.mpf(i[60:79])))
                    cnt += 1
                elif cnt == 2:
                    eph2.append((mp.mpf(i[3:22]),mp.mpf(i[22:41]),mp.mpf(i[41:60]),mp.mpf(i[60:79])))
                    cnt += 1
                elif cnt == 3:
                    eph3.append((mp.mpf(i[3:22]),mp.mpf(i[22:41]),mp.mpf(i[41:60]),mp.mpf(i[60:79])))
                    cnt += 1
                elif cnt == 4:
                    eph4.append((mp.mpf(i[3:22]),mp.mpf(i[22:41]),mp.mpf(i[41:60]),mp.mpf(i[60:79])))
                    cnt += 1
                elif cnt == 5:
                    eph5.append((mp.mpf(i[3:22]),mp.mpf(i[22:41]),mp.mpf(i[41:60]),mp.mpf(i[60:79])))
                    cnt += 1
                elif cnt == 6:
                    eph6.append((mp.mpf(i[3:22]),mp.mpf(i[22:41]),mp.mpf(i[41:60]),mp.mpf(i[60:79])))
                    cnt += 1
                elif cnt == 7:
                    eph7.append((mp.mpf(i[3:22]),mp.mpf(i[22:41]),mp.mpf(i[41:60]),mp.mpf(i[60:79])))
                    cnt = 1

    #--- name columns and put into Data Frames ---#
    eph = pd.DataFrame(eph,columns=('svid','year','month','day','hour','minutes','sec','af0','af1','af2'))
    eph1 = pd.DataFrame(eph1,columns=('IODE','crs','delta_n','M0'))
    eph2 = pd.DataFrame(eph2,columns=('Cuc','e','Cus','Asqrt'))
    eph3 = pd.DataFrame(eph3,columns=('Toe','Cic','OMEGA','Cis'))
    eph4 = pd.DataFrame(eph4,columns=('i0','Crc','omega','OMEGA_DOT'))
    eph5 = pd.DataFrame(eph5,columns=('IDOT','codeL2','GPS_Week','L2Pdata'))
    eph6 = pd.DataFrame(eph6,columns=('accuracy','health','TGD','IDOC'))
    eph7 = pd.DataFrame(eph7,columns=('ttx','Fit_interval','UNKNOWN1','UNKOWN2'))
    ## --- cat them all together, each line is 1 set of data ---- ## 
    ephf = pd.concat((eph,eph1,eph2,eph3,eph4,eph5,eph6,eph7),axis=1)
    
    return ephf
