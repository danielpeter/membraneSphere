#! /usr/bin/env python

# script for sorting a 'ttkernel.*' file to have increasing latitudes
# arguments:
#       $1 - filename e.g. 'Data/tmp/ttkernel.dat'

import os,math,sys
#import scipy

lon =[ ]
lat =[ ]
fileline=[ ]
commentline=[ ]

# get argument
if len(sys.argv) > 1:
  print '# argument 1:'+sys.argv[1] 
  filename=sys.argv[1]
else:
  print 'usage:'
  print '     sortLatitudes.py kernel.dat'
  print
  sys.exit(1)
  
# open data file
inputFile=filename
fileInput = open(inputFile, 'r')

# read file entries
commentline.append( '# reading file '+inputFile+'  ... ' )

index=0
line=fileInput.readline() 
while line :
  if len(line) < 10 :
    commentline.append( '# line too short:'+line )
    line=fileInput.readline()
    continue
  
  if line[1] == "#" :
    commentline.append( line[0:len(line)-1] )
    line=fileInput.readline()
    continue
     
  lineitems = line.split()

  x0=lineitems[0]
  x1=lineitems[1]
  x2=lineitems[2]

  lon.append( float(x0) )
  lat.append( float(x1) )
  fileline.append( line[0:len(line)-1] )
  index=index+1   

  line=fileInput.readline()

fileInput.close()
commentline.append( '# entries read in:'+str(index) )

# bubblesort
commentline.append( '# sorting entries ...' )
false=0
true=1
swapped= true
while( swapped ):
  swapped=false
  j=0
  while j < index-1:
    if lat[j+1]< lat[j]:
      swapped=true
      lat[j],lat[j+1]=lat[j+1],lat[j]
      fileline[j],fileline[j+1]=fileline[j+1],fileline[j]
    j=j+1
    
# print output    
for line in commentline:
  print line

for i in range(index):  
  print fileline[i]  

