# script for parameter search
import os
import scipy
import math

inputFile='simulation.0000.dat'
fileInput = open(inputFile, 'r')
count = 0
line=fileInput.readline() 

while line :
  if len(line) < 10 :
    print 'line too short'
    exit

  lineitems = line.split()

  x=lineitems[0]
  y=lineitems[1]
  z=lineitems[2]
  value=lineitems[3]

  # conversion
  x=float(x)
  y=float(y)
  z=float(z)
  value=float(value)

  longitude=math.atan2(y,x)
  latitude=math.asin(z)
  if longitude < 0 :
    longitude=math.pi+longitude

  longitude=longitude*180.0/math.pi
  latitude=latitude*180.0/math.pi


  print '%(lon) 12f %(lat) 12f %(val) 12f' %{"lon": longitude, "lat":latitude, "val":value} 

  line=fileInput.readline()
  count=count+1

fileInput.close()
 

