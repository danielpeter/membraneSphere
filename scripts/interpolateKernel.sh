#!/bin/bash

# ---------------------------------------------------------------------------------------------------------
# interpolates kernel values by gmt-routine
# ( gives readable kernel file '*.regular.xyz' for the inversion matrix construction - matrix.f in SUWA)
# ---------------------------------------------------------------------------------------------------------

# arguments:
#       $1 = data-file name, e.g. ttkernel.rot.dat

# check
if [ ! -s "$1" ]; then
  echo "usage: "
  echo "    interpolateKernel.sh kernel.dat"
  echo
  exit
fi

# ---------------------------------------------------------------------------------------------------------
# PARAMETERS
# ---------------------------------------------------------------------------------------------------------

region=-R0/360/-90/90
increment=-I1/1
latlon=no
crosssections=no

# input file
##datafilename=analytic.L150_90.avg0.9e-3.n10.dat
#datafilename=tt_kernel.cmplete.L150.dat
datafilename=$1

# output file
datafilename1=$1.regular.xyz

# ---------------------------------------------------------------------------------------------------------
# PRE-PROCESS
# ---------------------------------------------------------------------------------------------------------
# preprocess data file
if [ ! -s $datafilename ];then
     echo can not find data file $datafilename
     exit
fi

# extract longitude,latitude and kernelvalue to temporary file
echo '# data' > tmpData.dat

if [ "$latlon" = "yes" ];then
  gawk '{if( (substr($1,1,1)!="#") && ($1!="") && (substr($1,1,1)~"[0-9]") ) \
	print($2,$1,$3); }' $datafilename >> tmpData.dat
else
   gawk '{if( (substr($1,1,1)!="#") && ($1!="") && (substr($1,1,1)~"[0-9]") ) \
	print($1,$2,$3); }' $datafilename >> tmpData.dat
fi

# ---------------------------------------------------------------------------------------------------------
# GMT
# ---------------------------------------------------------------------------------------------------------

gmt gmtset HEADER_FONT_SIZE 14 MAP_ANNOT_OBLIQUE 0 BASEMAP_TYPE plain

# gmt interpolation
gmt blockmean tmpData.dat $region $increment > tmpData.med.xyz
gmt surface tmpData.med.xyz  $region $increment -Gimage.bin -T0.25 -C0.1

echo "#kernel values from" $datafilename > $datafilename1
echo "#lon lat val" >> $datafilename1
gmt grd2xyz image.bin >> $datafilename1

# cross-sections
if [ "$crosssections" = "yes" ];then
crosslon=45
echo '#cross-section at longitude ' $crosslon ' from' $datafilename > $datafilename.$crosslon.xyz
echo '#lon lat val' >> $datafilename.$crosslon.xyz
gawk '{if( (substr($1,1,1)!="#") && ($1!="") && (substr($1,1,1)~"[0-9]") ) \
	if( $1=='$crosslon') print($1,$2,$3); }' $datafilename1 >> $datafilename.$crosslon.xyz

crosslon=80
echo '#cross-section at longitude ' $crosslon ' from' $datafilename > $datafilename.$crosslon.xyz
echo '#lon lat val' >> $datafilename.$crosslon.xyz
gawk '{if( (substr($1,1,1)!="#") && ($1!="") && (substr($1,1,1)~"[0-9]") ) \
	if( $1=='$crosslon') print($1,$2,$3); }' $datafilename1 >> $datafilename.$crosslon.xyz

crosslon=70
echo '#cross-section at longitude ' $crosslon ' from' $datafilename > $datafilename.$crosslon.xyz
echo '#lon lat val' >> $datafilename.$crosslon.xyz
gawk '{if( (substr($1,1,1)!="#") && ($1!="") && (substr($1,1,1)~"[0-9]") ) \
	if( $1=='$crosslon') print($1,$2,$3); }' $datafilename1 >> $datafilename.$crosslon.xyz

crosslon=20
echo '#cross-section at longitude ' $crosslon ' from' $datafilename > $datafilename.$crosslon.xyz
echo '#lon lat val' >> $datafilename.$crosslon.xyz
gawk '{if( (substr($1,1,1)!="#") && ($1!="") && (substr($1,1,1)~"[0-9]") ) \
	if( $1=='$crosslon') print($1,$2,$3); }' $datafilename1 >> $datafilename.$crosslon.xyz

crosslon=10
echo '#cross-section at longitude ' $crosslon ' from' $datafilename > $datafilename.$crosslon.xyz
echo '#lon lat val' >> $datafilename.$crosslon.xyz
gawk '{if( (substr($1,1,1)!="#") && ($1!="") && (substr($1,1,1)~"[0-9]") ) \
	if( $1=='$crosslon') print($1,$2,$3); }' $datafilename1 >> $datafilename.$crosslon.xyz

echo "cross-section at longitude 10/20/45/70/80 to $datafilename.--.xyz"
fi

echo "************************************************************************************"
echo "grid plotted to file: $datafilename1"
linenumbers=`gawk 'END { print NR }' $datafilename1`
echo "number of lines = $linenumbers ( 2 comment lines + 181*361 = 65343 )"
if [ $linenumbers != 65343 ]; then
  echo "uncorrect kernel!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  pause
  exit
fi
echo "************************************************************************************"
echo



