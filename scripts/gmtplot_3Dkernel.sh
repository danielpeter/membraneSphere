#!/bin/bash

# ---------------------------------------------------------------------------------------------------------
# gmtplot_3Dkernel.sh
#
# draws the kernel values from a given data file on a gmt map
# ---------------------------------------------------------------------------------------------------------

datafilename=$1

# check
if [ ! -s "$1" ]; then
  echo "usage: "
  echo "    gmtplot_3Dkernel.sh kernel.dat"
  echo
  exit
fi

# ---------------------------------------------------------------------------------------------------------
# PARAMETERS
# ---------------------------------------------------------------------------------------------------------

#datafilename=analytic.L150_90.avg0.9e-3.n10.dat
#zdepth=-Jz0.745
#datafilename=tt_kernel.cmplete.L150.dat

# parameters
projection=-JQ0/15
Plotregion=-R-10/100/-60/60
region=-R0/360/-90/90
offsets=  #'-X1.5 -Y14'
portrait=-P
colormap=seis.L150.4.cpt
coastoffset=-Y8  #-Y13.85
zdepth=-Jz0.22   #-Jz1.0 #-Jz0.732
perspective=-E175/35

#title='travel time sensitivity kernel'

# check if colormap file exists in ../../scripts/ folder
if [ -e ../../scripts/$colormap ]; then colormap=../../scripts/$colormap; fi
if [ ! -e $colormap ]; then
  echo ""; echo "Error: color map not found: $colormap"; echo "Please find colormaps in scripts/ folder"; echo "";
  exit 1
fi
echo "using color map: " $colormap


ps_filename=$datafilename.3D.ps
jpg_filename=$datafilename.3D.jpg

echo
echo "plotting to file: " $ps_filename
echo

# ---------------------------------------------------------------------------------------------------------
# PRE-PROCESS
# ---------------------------------------------------------------------------------------------------------

# check if polygon data available
if [ ! -s $datafilename ];then
  echo "can not find polygon data file $datafilename"
  exit 1
fi

# ---------------------------------------------------------------------------------------------------------
# GMT
# ---------------------------------------------------------------------------------------------------------

gmt gmtset PS_MEDIA letter HEADER_FONT_SIZE 14 MAP_ANNOT_OBLIQUE 0 BASEMAP_TYPE plain GMT_HISTORY false

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

gmt xyz2grd $datafilename -Gimage.bin -I1/1 $region  -N0.0 -V

#gmt grdimage image.bin $portrait $offset $region -G $projection -Bg15/g15 -T -K -V -C$colormap > $ps_filename

gmt grdview image.bin $perspective $zdepth -Qs $portrait $Plotregion $projection -K -V -C$colormap > $ps_filename
#gmt psscale -D15/3.5/4/0.5 -Ba.2f0.1:"": -C$colormap -O -Y-1.5 $portrait -V -K  >> $ps_filename
#gmt pscoast $coastoffset $Plotregion $projection $portrait $perspective -W1 -Dc -A10000 -O -K >> $ps_filename
gmt psbasemap $coastoffset $Plotregion $verbose $projection $portrait $perspective -Bg15/g15:."": -O   >> $ps_filename

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo "plotted to $ps_filename"


# to play with...

# Use pscoast to plot a map with different colors for land and sea rivers(I1), political(N1)
#gmt pscoast    $region $projection -Di -W1/255/255/255 -A1000  \
#         -S50/50/50 -G0/0/0 -N1/255/255/255 $portrait \
#         -K   $offsets > $ps_filename

#gmt xyz2grd timelag.L150.l6.left.dat -Gimage.bin -I1/1 -R-50/-6/-90/90  -N0.0 -V
#gmt grdimage image.bin $portrait $offset $region -G $projection -Bg15/g15 -T -K -V -C$colormap -K -O >> $ps_filename
# Add custom xy-data from $datafilename
# take size 0.4/0.05 for grid level 4/8 and Europe
# take size 0.1      for grid level 4 and World
#gawk '{if((substr($1,1,1)!="#")&&($1!=""))print($1,$2,$3,0.08)}' $datafilename | \
#          $gmtbin/psxy $verbose $region $projection -O -K -G002/170/002 -Sc -C$colormap >> $ps_filename
#gawk '{if((substr($1,1,1)!="#")&&($1!=""))print((90-$1),$2,$3,0.08)}' $datafilename | \
#          $gmtbin/psxy $verbose $region $projection -O -K -G002/170/002 -Sc -C$colormap >> $ps_filename
#gawk '{if((substr($1,1,1)!="#")&&($1!=""))print((90-$1),-$2,$3,0.08)}' $datafilename | \
#          $gmtbin/psxy $verbose $region $projection -O -K -G002/170/002 -Sc -C$colormap >> $ps_filename
#gawk '{if((substr($1,1,1)!="#")&&($1!=""))print($1,-$2,$3,0.08)}' $datafilename | \
#          $gmtbin/psxy $verbose $region $projection -O -K -G002/170/002 -Sc -C$colormap >> $ps_filename
## color scale legend -B3:"sensitivity":
#gmt psscale  $portrait -C$colormap -D15/3.5/4/.5 -O -L -B20:"":  >> $ps_filename
# Use psbasemap for basic map layout, possible title
#     and complete the plot    #


echo
echo "converting image..."

# pdf
#convert $ps_filename $pdf_filename
#if [ -s $pdf_filename ]; then
#  rm -f $ps_filename
#fi
#open $pdf_filename

# jpg
magick -density 300 $ps_filename -quality 90 $jpg_filename

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo "*******************************************************************"
echo
echo "image plotted to file           : $jpg_filename"
echo
echo "*******************************************************************"
