#!/bin/bash


# check
if [ ! -s "$1" ]; then
  echo "usage: "
  echo "    gmtplot_3Dkernel.sh kernel.dat"
  echo
  exit
fi

#datafilename=analytic.L150_90.avg0.9e-3.n10.dat
#zdepth=-Jz0.745
#datafilename=tt_kernel.cmplete.L150.dat

datafilename=$1
ps_filename=$1.3D.ps


# parameters
projection=-JQ0/15
Plotregion=-R-10/100/-60/60
region=-R0/360/-90/90
offsets=  #'-X1.5 -Y14'
portrait=-P
colormap=seis.L150.4.cpt
coastoffset=-Y13.45  #13.85
zdepth=-Jz1.0       #0.732
perspective=-E175/35



title='travel time anomaly kernel - Love 150 s'


# check if polygon data available
if [ ! -s $datafilename ];then
     echo "can not find polygon data file $datafilename"
     exit
fi

gmt gmtset HEADER_FONT_SIZE 14 MAP_ANNOT_OBLIQUE 0 BASEMAP_TYPE plain GMT_HISTORY false
gmt xyz2grd $datafilename -Gimage.bin -I1/1 $region  -N0.0 -V

#gmt grdimage image.bin $portrait $offset $region -G $projection -Bg15/g15 -T -K -V -C$colormap > $ps_filename

gmt grdview image.bin $perspective $zdepth -Qs -P $Plotregion $projection -K -V -C$colormap > $ps_filename
gmt psscale -D15/3.5/4/0.5 -Ba.2f0.1:"": -C$colormap -O -Y-1.5 $portrait -V -K  >> $ps_filename
gmt pscoast $coastoffset $Plotregion $projection $portrait $perspective -W1 -Dc -A10000 -O -K >> $ps_filename
gmt psbasemap  $Plotregion $verbose $projection $perspective -Bg15/g15:."": -O   >> $ps_filename

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

