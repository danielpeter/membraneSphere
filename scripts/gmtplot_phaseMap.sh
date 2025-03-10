#!/bin/bash

# arguments:
#  $1 = phase velocity map , e.g. PhaseMap.dat
#  $2 = show relative percentages
#  $3 = reference velocity

# check
if [ ! -s "$1" ]; then
  echo "usage: "
  echo "    gmtplot_phaseMap PhaseMap.dat [yes] [c0]"
  echo
  exit
fi

# input arguments
datafilename=$1
# relative percent of a reference speed
percent=$2
# reference speed
c0=$3

## GMT
# Geographical variables:
projection=-JQ0/18
#projection=-JM15

#region=-R-180/180/-70/70
region=-R0/360/-90/90

offsets='-X1.5 -Y14'
portrait=-P
verbose=-V

# colormap / annotation
if [ "$percent" = "yes" ];then
  # colortable
  colormap=seis.L150.1.cpt
  scaleAnotate=-Ba1:"%":
else
  # colortable
  colormap=seis.L150.3.cpt
  scaleAnotate=-Ba.1f1:"km/s":
fi

# relative percent of reference speed
c_ref=`grep 'reference velocity' $datafilename | sed -E 's/.*=\s*([0-9.]+).*/\1/'`
if [ "${c_ref}" != "" ]; then
  echo "data file  : $datafilename"
  echo "             file header has reference velocity = ${c_ref}"
  echo
  # check if phase map is giving perturbations in percent
  if [[ "$datafilename" == *".percent."* ]] && [[ "${c_ref}" != "" ]]; then
    echo "data file is phase map given in percent."
    # colortable
    colormap=seis.L150.1.cpt
    scaleAnotate=-Ba1:"%":
    echo
    # plotting as perturbation to reference velocity in percent only makes sense for phase maps giving absolute phase speed values
    if [ "$percent" == "yes" ]; then
      echo "Phase map is already given in percent - calculating perturbations to reference velocity in percent not needed."
      exit 1
    fi
  fi
fi

if [ "$percent" = "yes" ];then
  # relative percent of reference speed
  echo
  echo "plotting phase speed perturbation as percentages"
  if [ "${c0}" == "" ]; then c0=4.77915; fi
  echo "         relative speed percentages of c0=${c0} ..."
  echo
fi


# check if colormap file exists in ../../scripts/ folder
if [ -e ../../scripts/$colormap ]; then colormap=../../scripts/$colormap; fi
if [ ! -e $colormap ]; then
  echo ""; echo "Error: color map not found: $colormap"; echo "Please find colormaps in scripts/ folder"; echo "";
  exit 1
fi

# file output
ps_filename=$1.ps
title='phase velocity map'

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
     echo "can not find polygon data file $datafilename"
     exit
fi

# Add custom xy-data from $datafilename
# take size 0.4/0.05 for grid level 4/8 and Europe
# take size 0.1      for grid level 4 and World

# create temporary table
if [ "$percent" = "yes" ];then
  gawk '{if((substr($1,1,1)!="#")&&($1!=""))print($1,$2,(($3-c0)/c0*100),0.05)}' c0=${c0} $datafilename > gmt.table.txt
else
  gawk '{if((substr($1,1,1)!="#")&&($1!=""))print($1,$2,$3,0.05)}' $datafilename > gmt.table.txt
fi

# infos
echo "extracted table infos:"
gmt info gmt.table.txt
echo

# plotting
gmt psxy gmt.table.txt $verbose $region $projection -K -G002/170/002 -Sc -C$colormap $offsets $portrait > $ps_filename

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
gmt pscoast $region $projection -Di -A1000 -W1 -O -K $verbose >> $ps_filename

## color scale legend -B3:"speed":
gmt psscale -D9/-1/4/0.5h $scaleAnotate -C$colormap -V -O -K >> $ps_filename

# Use psbasemap for basic map layout, possible title
#     and complete the plot
#$gmtbin/psbasemap  $region $verbose $projection \
#	  -Ba30f10.000/a30f10.000:."$title":WeSn \
#	  -O   >> $ps_filename
gmt psbasemap $projection $region  -Ba90f90  -O >> $ps_filename

echo converting...
# jpg
magick -density 300 $1.ps -quality 90 $1.jpg
#rm -f $1.ps

echo
echo "image plotted to: $1.jpg"
echo


