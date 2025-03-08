#!/bin/bash

# arguments:
#  $1 = phase velocity map , e.g. PhaseMap.dat
#  $2 = show relative percentages

# check
if [ ! -s "$1" ]; then
  echo "usage: "
  echo "    gmtplot_phaseMap PhaseMap.dat [yes]"
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

# relative percent of a reference speed
percent=$2

if [ "$percent" = "yes" ];then
  # relative percent of reference speed
  echo "plotting relative speed percentages of c0=4.77915 ..."

  # colortable
  colormap=seis.L150.1.cpt

  scaleAnotate=-Ba1:"%":
else
  # colortable
  colormap=seis.L150.3.cpt

  scaleAnotate=-Ba.1f1:"km/s":
fi


# file output
datafilename=$1
ps_filename=$1.ps
title='phase velocity map'

################################################################################
#
# GMT plotting commands follow
#
################################################################################

gmtset HEADER_FONT_SIZE 14 OBLIQUE_ANOTATION 0 BASEMAP_TYPE plain

#gmtset GLOBAL_X_SCALE 1. GLOBAL_Y_SCALE 1.
#gmtset ANOT_FONT Times-Roman ANOT_FONT_SIZE 15
#gmtset LABEL_FONT Times-Roman LABEL_FONT_SIZE 18
#gmtset HEADER_FONT Times-Bold HEADER_FONT_SIZE 18
#gmtset DEGREE_FORMAT 1 DOTS_PR_INCH 300 MEASURE_UNIT cm
#gmtset BASEMAP_TYPE fancy



# check if polygon data available
if [ ! -s $datafilename ];then
     echo can not find polygon data file $datafilename
     exit
fi

# Add custom xy-data from $datafilename
# take size 0.4/0.05 for grid level 4/8 and Europe
# take size 0.1      for grid level 4 and World


if [ "$percent" = "yes" ];then
  gawk '{if((substr($1,1,1)!="#")&&($1!=""))print($1,$2,(($3-4.77915)*100/4.77915),0.05)}' $datafilename | \
	 psxy $verbose $region $projection -K -G002/170/002 -Sc -C$colormap $offsets $portrait > $ps_filename
else
  gawk '{if((substr($1,1,1)!="#")&&($1!=""))print($1,$2,$3,0.05)}' $datafilename | \
	 psxy $verbose $region $projection -K -G002/170/002 -Sc -C$colormap $offsets $portrait > $ps_filename
fi


# source
  #gawk '{if((substr($1,1,1)!="#")&&($1!=""))if(NR==77607)print($1,$2,0.5)}' $datafilename | \
  #	  psxy $verbose $region $projection -O -K -W0/0/0 -G20/20/20 -Sa >> $ps_filename #-C$colormap
  #echo source:
  #gawk '{if((substr($1,1,1)!="#")&&($1!=""))if(NR==77607)print($1,$2,0.5)}' $datafilename > coords

# receiver
  #gawk '{if((substr($1,1,1)!="#")&&($1!=""))if(NR==59)print($1,$2,0.5)}' $datafilename | \
  #	  psxy $verbose $region $projection -O -K -W0/0/0 -G20/20/20 -St >> $ps_filename
  #echo receiver:
  #gawk '{if((substr($1,1,1)!="#")&&($1!=""))if(NR==59)print($1,$2,0.5)}' $datafilename >> coords
  #psxy -W5/220/220/220  $projection $region coords -V -O -K >> $ps_filename

# Use pscoast to plot a map with different colors for land and sea rivers(I1), political(N1)
#  -W1/255/255/255 -S50/50/50 -N1/255/255/255 -G0/0/0
pscoast $region $projection -Di -A1000 -W1 -O -K $verbose >> $ps_filename

## color scale legend -B3:"speed":
psscale -D9/-1/4/0.5h $scaleAnotate -C$colormap -V -O -K >> $ps_filename

# Use psbasemap for basic map layout, possible title
#     and complete the plot
#$gmtbin/psbasemap  $region $verbose $projection \
#	  -Ba30f10.000/a30f10.000:."$title":WeSn \
#	  -O   >> $ps_filename
psbasemap $projection $region  -Ba90f90  -O >> $ps_filename

echo converting...
convert $1.ps $1.pdf

echo
echo image plotted to: $1.pdf
echo


