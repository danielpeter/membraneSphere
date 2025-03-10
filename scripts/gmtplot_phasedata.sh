#!/bin/bash

# arguments:
#  $1 = grid phase velocities e.g. phase.6.L150.dat, L150.gsh-012.xyz

# check
if [ ! -s "$1" ]; then
  echo "usage: "
  echo "    gmtplot_phasedata L150.gsh-012.xyz"
  echo
  exit
fi

# Geographical variables:
projection=-JQ0/18
#projection=-JM15

#region=-R-180/180/-70/70
region=-R0/360/-90/90

offsets='-X1.5 -Y14'
portrait=-P
verbose=-V

colormap=seis.L150.2.cpt
datafilename=$1


title='phase velocity map'
ps_filename=tmp.ps

################################################################################
#
# GMT plotting commands follow
#
################################################################################

gmt gmtset HEADER_FONT_SIZE 14 MAP_ANNOT_OBLIQUE 0 BASEMAP_TYPE plain GMT_HISTORY false

#gmt gmtset GLOBAL_X_SCALE 1. GLOBAL_Y_SCALE 1.
#gmt gmtset ANOT_FONT Times-Roman ANOT_FONT_SIZE 15
#gmt gmtset LABEL_FONT Times-Roman LABEL_FONT_SIZE 18
#gmt gmtset HEADER_FONT Times-Bold HEADER_FONT_SIZE 18
#gmt gmtset DEGREE_FORMAT 1 DOTS_PR_INCH 300 MEASURE_UNIT cm
#gmt gmtset BASEMAP_TYPE fancy



# check if polygon data available
if [ ! -s $datafilename ];then
     echo can not find polygon data file $datafilename
     exit
fi

# Add custom xy-data from $datafilename
# take size 0.4/0.05 for grid level 4/8 and Europe
# take size 0.1      for grid level 4 and World
gawk '{if((substr($1,1,1)!="#")&&($1!=""))print($1,$2,$3,0.05)}' $datafilename | \
	 gmt psxy $verbose $region $projection -K -G002/170/002 -Sc -C$colormap $offsets $portrait > $ps_filename

# source
gawk '{if((substr($1,1,1)!="#")&&($1!=""))if(NR==77607)print($1,$2,0.5)}' $datafilename | \
	  gmt psxy $verbose $region $projection -O -K -W0/0/0 -G20/20/20 -Sa >> $ps_filename #-C$colormap
#echo source:
gawk '{if((substr($1,1,1)!="#")&&($1!=""))if(NR==77607)print($1,$2,0.5)}' $datafilename > coords

# receiver
gawk '{if((substr($1,1,1)!="#")&&($1!=""))if(NR==59)print($1,$2,0.5)}' $datafilename | \
	  gmt psxy $verbose $region $projection -O -K -W0/0/0 -G20/20/20 -St >> $ps_filename
#echo receiver:
gawk '{if((substr($1,1,1)!="#")&&($1!=""))if(NR==59)print($1,$2,0.5)}' $datafilename >> coords

gmt psxy -W5/220/220/220  $projection $region coords -V -O -K >> $ps_filename

# Use pscoast to plot a map with different colors for land and sea rivers(I1), political(N1)
#  -W1/255/255/255 -S50/50/50 -N1/255/255/255 -G0/0/0
gmt pscoast $region $projection -Di -A1000 -W1 -O -K  >> $ps_filename

## color scale legend -B3:"speed":
gmt psscale -D9/-1/4/0.5h -Ba.1f1:"km/s": -C$colormap -V -O -K >> $ps_filename

# Use psbasemap for basic map layout, possible title
#     and complete the plot
#$gmtbin/psbasemap  $region $verbose $projection \
#	  -Ba30f10.000/a30f10.000:."$title":WeSn \
#	  -O   >> $ps_filename
gmt psbasemap $projection $region  -Ba90f90  -O >> $ps_filename

echo converting...

# pdf
#convert tmp.ps tmp.pdf

# jpg
magick -density 300 $tmp.ps -quality 90 tmp.jpg

rm -f tmp.ps

echo
echo "image plotted to tmp.jpg"
echo
