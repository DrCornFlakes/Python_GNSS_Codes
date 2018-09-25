# -*- coding: utf-8 -*-
"""
Created on Mon Jul 02 10:44:27 2018

@author: Dawn Merriman

"""

import pandas as pd
import csv 

def SimpleLineGrep(ASCIIfileName,filterKeyword):
  from cStringIO import StringIO
  s = StringIO()

  with open(ASCIIfileName) as f:
      for line in f:
          if line.startswith(filterKeyword):
              s.write(line)
  s.seek(0) # "rewind" to the beginning of the StringIO object
  # slow, test dump implementation
  # outFile = open("%s/dmp.csv" % ASCIIfileName[:ASCIIfileName.index('/')], "w")
  # outFile.write(s.getvalue())
  # outFile.close()

  return s

# reads GNSSLogger log file into pandas data frame using Google Logger 2.0
def readGNSSLogger (data_file):
  print("filtering PR from %s" % data_file)
  RawMeas = SimpleLineGrep(data_file,'Raw')
  # for 7.0 colNames = ["Raw","ElapsedRealtimeMillis","TimeNanos","LeapSecond","TimeUncertaintyNanos","FullBiasNanos","BiasNanos","BiasUncertaintyNanos","DriftNanosPerSecond","DriftUncertaintyNanosPerSecond","HardwareClockDiscontinuityCount","Svid","TimeOffsetNanos","State","ReceivedSvTimeNanos","ReceivedSvTimeUncertaintyNanos","Cn0DbHz","PseudorangeRateMetersPerSecond","PseudorangeRateUncertaintyMetersPerSecond","AccumulatedDeltaRangeState","AccumulatedDeltaRangeMeters","AccumulatedDeltaRangeUncertaintyMeters","CarrierFrequencyHz","CarrierCycles","CarrierPhase","CarrierPhaseUncertainty","MultipathIndicator","SnrInDb","ConstellationType","AgcDb","CarrierFrequencyHz"]
  # for 8.0
  colNames = ["Raw","ElapsedRealtimeMillis","TimeNanos","LeapSecond","TimeUncertaintyNanos","FullBiasNanos","BiasNanos","BiasUncertaintyNanos","DriftNanosPerSecond","DriftUncertaintyNanosPerSecond","HardwareClockDiscontinuityCount","Svid","TimeOffsetNanos","State","ReceivedSvTimeNanos","ReceivedSvTimeUncertaintyNanos","Cn0DbHz","PseudorangeRateMetersPerSecond","PseudorangeRateUncertaintyMetersPerSecond","AccumulatedDeltaRangeState","AccumulatedDeltaRangeMeters","AccumulatedDeltaRangeUncertaintyMeters","CarrierFrequencyHz","CarrierCycles","CarrierPhase","CarrierPhaseUncertainty","MultipathIndicator","SnrInDb","ConstellationType","AgcDb","CarrierFrequencyHz"]
  dataFrame = pd.read_csv(RawMeas, delimiter = ",",error_bad_lines=False,header=None,
                          usecols=range(1,len(colNames)),names= colNames,
                          encoding = 'utf-8-sig',na_values = ["NULL",""],engine ='c')
  return dataFrame

def readgnss(data_file):
    d = {}
    time = []
    with open(data_file, 'rb') as csvfile:
        unpack = csv.reader(csvfile, delimiter=",", quotechar='|')
        for row in unpack:
            if "Raw" in row:
                if float(row[18]) > abs(2*float(row[17])):
                    pass  #if uncertainity is double the value
                else:
                    d.setdefault((row[28],row[11]),[]).append((row[1],row[17],row[18]))
                    time.append(float(row[1]))

    min_time = min(time)
    
    return min_time,d

def readgnss2(data_file):
    d = {}
    time = []
    with open(data_file, 'rb') as csvfile:
        unpack = csv.reader(csvfile, delimiter=",", quotechar='|')
        for row in unpack:
            if "Raw" in row:
                if float(row[18]) > abs(2*float(row[17])):
                    pass  #if uncertainity is double the value
                else:
                    d.setdefault((row[28],row[11]),[]).append((row[1],row[17],row[18]))
                    time.append(float(row[1]))

    min_time = min(time)
    
    return min_time,d