#!/usr/bin/env python
#
# script for parameter search
#
# varies either the delta location or the perturbation size in the 'Parameter_Input' file
# and calls the executable 'phaseshift'
#
# outputs the results on console

import os
import sys
import fileinput
import string

# make a backup of parameter file
inputFile='../Parameter_Input.org'
newFile='../Parameter_Input'
os.system('cp ../Parameter_Input ../Parameter_Input.org')

# change delta location
DELTALOCATION=0

# change delta perturbation size
STRONGNESS=1


# function: replaceParameter 
# the next line of the given file after the given parametername is changed accordingly
def replaceParameter(inputFile,newFile,parametername,newvalue):
    os.system('cp '+inputFile+' tmp.dat')
    fileToChange = open(newFile, 'w')
    linenumber = 1
    nextline = -1
    for line in fileinput.input( 'tmp.dat' ):
      #if string.find(line, parametername) >= 0:
      #  nextline = linenumber+1
      #if nextline == linenumber:
      #  fileToChange.write(str(newvalue)+'\n')
      #  nextline = 0
      #  #print "replaced:",parametername,newvalue
      if string.find(line, parametername) >= 0:
        nextline = linenumber
        length = len(line)
        #print "line length:",length
        #print "    ",line
        newline = line[0:34] + " " + str(newvalue) + "\n"
        #print "    ",newline
        fileToChange.write(newline)
      else:
        fileToChange.write(line)
      linenumber = linenumber + 1
      
    # close files      
    fileinput.close()
    fileToChange.close()  
    os.system('rm tmp.dat')
    
    # check if new parameter is set
    if nextline == -1:
      print "COULD NOT REPLACE PARAMETER VALUE:",parametername,newvalue
    return


# change delta location
if DELTALOCATION == 1:
  #starting values
  deltaLatitude=0
  deltaLongitude=45
  #increment
  dLat = 1
  # loop latitude
  while deltaLatitude <= 90:
    print 'deltaLocation:'+str(deltaLatitude)+' '+str(deltaLongitude)
    # create new input parameter file
    fileToChange = open(newFile, 'w')
    # replace the Delta location
    linenumber = 1
    nextline = -1
    for line in fileinput.input( inputFile ):
      if line == 'DLOCATION = \n':
        nextline = linenumber+1
      if nextline == linenumber:
        fileToChange.write(str(deltaLatitude)+' '+str(deltaLongitude)+'\n')
        nextline = -1
      else:
        fileToChange.write(line)
      linenumber = linenumber + 1    
    fileToChange.close()  
    
    #print 'calling: ./propagation > Data/delta.dr8/delta.dr8strong.narrow.'+str(deltaLatitude)+'.dat'
    print 'calling: mpirun -np 2 propagation > seismo.L300.'+str(deltaLatitude)+'.dat'
    os.system('mpirun -np 2 propagation > seismo.L300.'+str(deltaLatitude)+'.dat') 
    deltaLatitude = deltaLatitude+dLat
    
    
if STRONGNESS == 1:
  #### perturbation size in percent
  # long test
  #percents=[-20,-10,-8,-6,-4,-2,-1,-0.5,-0.3,-0.1,-0.05,-0.01,-0.001,0.0,0.001,0.01,0.05,0.1,0.3,0.5,1,2,4,6,8,10,20]
  # short test
  percents=[-20,-10,-5,-2,-0.2,0.0,0.2,2,5,10,20]
  
  #### simulation parameters
  # love wave phase velocity
  phasevelocityRef=4.77915   # L300: 5.22906, L600: 6.36841
  # change 'Parameter_Input' file accordingly 
  replaceParameter(inputFile,newFile,'CPHASE','L150')
  replaceParameter(newFile,newFile,'DELTA','.true.')
    
  # only single location
  replaceParameter(newFile,newFile,'MOVEDELTA','.false.')
  # scatterer location
  replaceParameter(newFile,newFile,'DLOCATION','0 45')
  # other parameters
  replaceParameter(newFile,newFile,'PARALLELSEISMO','.true.')
  replaceParameter(newFile,newFile,'SECONDDELTA','.false.')
  replaceParameter(newFile,newFile,'MANYRECEIVERS','.false.')
  replaceParameter(newFile,newFile,'MANYKERNELS','.false.')
  
  #### output
  print '# phase velocity perturbation effects'
  print '# lat  lon kernel  relativePhaseshift  perturbationPercent phasevelocityIncrement  phasevelocityofPerturbation'
  for perc in percents:
    #### calculate perturbation in [km/s]
    phasevelocity=phasevelocityRef+phasevelocityRef*perc/100.0
    phasevelocityIncrement=phasevelocity-phasevelocityRef
    print "# perturbation in percent:",perc
    print "# perturbation in km/s:",phasevelocityIncrement
    print "# perturbations phase velocity:",phasevelocity
    
    #### new perturbation size
    replaceParameter(newFile,newFile,'DPERTURBATION',phasevelocityIncrement)
    
    #debug
    #os.system('cat Parameter_Input')
    
    #print 'calling : mpirun -np 1 phaseshift > phaseshift.perc.'+str(perc)+'.log'
    os.system('mpirun -np 1 phaseshift > phaseshift.perc.'+str(perc)+'.log') 
    
    #### print phaseshift value
    phaseshift=0.0
    for line in fileinput.input('OUTPUT/ttkernel.dat'):
      #print line
      if len(line) <= 1:
        continue
      if line[1] != '#':
        numbers=string.split(line)
        #latitude = float(numbers[1])
        #longitude= float(numbers[2])
        #phaseshift=float(numbers[3])
        latitude = numbers[0]
        longitude= numbers[1]
        kernel=numbers[2]
        phaseshift=numbers[3]
        #latitude = '%(#) 12e' % {"#":timelag}
        percentage= '%(#) 10.3f' %{"#":perc}
        phasevelocityInc= '%(#) 10.3f' %{"#":phasevelocityIncrement}
        
        # output
        print latitude,longitude,kernel,phaseshift,percentage,phasevelocityInc,phasevelocity
    
# backup original            
os.system('cp Parameter_Input.org Parameter_Input')            
