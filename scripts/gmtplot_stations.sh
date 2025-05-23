#!/bin/bash

# arguments:
#  $1 = station file , e.g. tmpReceiverStations001.dat
datafilename=$1

# check
if [ ! -s "$1" ]; then
  echo "usage: "
  echo "    gmtplot_stations tmpReceiverStations001.dat"
  echo
fi

# Geographical variables:
projection=-JQ0/18
#projection=-JM15

#region=-R-180/180/-70/70
region=-R0/360/-90/90

offsets='-X1.5 -Y14'
portrait=-P
verbose=-V

title='phase velocity map'
ps_filename=$datafilename.ps

################################################################################
#
# GMT plotting commands follow
#
################################################################################

gmt gmtset PS_MEDIA letter HEADER_FONT_SIZE 14 MAP_ANNOT_OBLIQUE 0 BASEMAP_TYPE plain GMT_HISTORY false

#gmt gmtset GLOBAL_X_SCALE 1. GLOBAL_Y_SCALE 1.
#gmt gmtset ANOT_FONT Times-Roman ANOT_FONT_SIZE 15
#gmt gmtset LABEL_FONT Times-Roman LABEL_FONT_SIZE 18
#gmt gmtset HEADER_FONT Times-Bold HEADER_FONT_SIZE 18
#gmt gmtset DEGREE_FORMAT 1 DOTS_PR_INCH 300 MEASURE_UNIT cm
#gmt gmtset BASEMAP_TYPE fancy

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# check if polygon data available
if [ ! -s $datafilename ];then
     echo can not find data file $datafilename
     exit
fi

# Add custom xy-data from $datafilename
# take size 0.4/0.05 for grid level 4/8 and Europe
# take size 0.1      for grid level 4 and World
gawk '{if((substr($1,1,1)!="#")&&($1!=""))print($1,$2,0.05)}' $datafilename | \
	 gmt psxy $verbose $region $projection -K -G002/170/002 -Sc $offsets $portrait -: > $ps_filename

# source
gawk '{if((substr($1,1,1)!="#")&&($1!=""))if(NR==1)print($1,$2,0.5)}' $datafilename | \
	  gmt psxy $verbose $region $projection -O -K -W0/0/0 -G20/20/20 -Sa -: >> $ps_filename #-C$colormap
#echo source:
gawk '{if((substr($1,1,1)!="#")&&($1!=""))if(NR==1)print($1,$2,0.5)}' $datafilename > coords

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

gmt psxy -W5/220/220/220  $projection $region coords -V -O -K -: >> $ps_filename

# Use pscoast to plot a map with different colors for land and sea rivers(I1), political(N1)
#  -W1/255/255/255 -S50/50/50 -N1/255/255/255 -G0/0/0
gmt pscoast $region $projection -Di -A1000 -W1 -O -K  >> $ps_filename

## color scale legend -B3:"speed":
#gmt psscale -D9/-1/4/0.5h -Ba.1f1:"km/s": -C$colormap -V -O -K >> $ps_filename

# Use psbasemap for basic map layout, possible title
#     and complete the plot
#gmt psbasemap  $region $verbose $projection \
#	  -Ba30f10.000/a30f10.000:."$title":WeSn \
#	  -O   >> $ps_filename
gmt psbasemap $projection $region  -Ba90f90  -O >> $ps_filename

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# pdf
#convert $ps_filename $datafilename.pdf
#open $datafilename.pdf

# jpg
magick -density 300 $ps_filename -quality 90 $datafilename.jpg

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

rm -f $ps_filename

echo
echo "image plotted to $datafilename.jpg"
echo

