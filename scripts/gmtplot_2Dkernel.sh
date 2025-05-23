#!/bin/bash

# ---------------------------------------------------------------------------------------------------------
# gmtplot_2Dkernel.sh
#
# draws the kernel values from a given data file on a gmt map
# ---------------------------------------------------------------------------------------------------------

# arguments:
#       $1 = data-file name, e.g. ttkernel.rot.dat
datafilename=$1
diff=$2

# check
if [ ! -s "$1" ]; then
  echo "usage: "
  echo "    gmtplot_2Dkernel.sh kernel.dat [diff]"
  echo
  exit
fi

# ---------------------------------------------------------------------------------------------------------
# PARAMETERS
# ---------------------------------------------------------------------------------------------------------

# interpolate grid data by gmt's spherical harmonics
interpolate=yes

# close-ups
#Plotregion=-R-10/100/-60/60
#Plotregion=-R-10/150/-60/60
#Plotregion=-R-60/150/-90/90
#Plotregion=-R35/175/-10/85 # tibet kernels

# full region
Plotregion=-R-180/180/-90/90

# map
projection=-JQ0/15
region=-R0/360/-90/90

# orientation
offsets='-X1.5 -Y14'
portrait=-P
latlon=no
scalePosition=-D7.5/0/6/0.5h
#scalePosition=-D5/0/6/0.5h
#scalePosition=-D12/0/6/0.5h # tibet kernels

#zdepth=-Jz2
#zdepth=-Jz0.745
#coastoffset=-Y13.85
#scaleoffset=-Y-1.5

##datafilename=analytic.L150_90.avg0.9e-3.n10.dat
#datafilename=tt_kernel.cmplete.L150.dat

if [ "$diff" = "diff" ];then
  colormap=seis.L150.diff1.cpt
  scaleAnotate=-Ba1:"":
else
  colormap=seis.L150.4.cpt
  scaleAnotate=-Ba5:"":
fi

# check if colormap file exists in ../../scripts/ folder
if [ -e ../../scripts/$colormap ]; then colormap=../../scripts/$colormap; fi
if [ ! -e $colormap ]; then
  echo ""; echo "Error: color map not found: $colormap"; echo "Please find colormaps in scripts/ folder"; echo "";
  exit 1
fi
echo "using color map: " $colormap

#title='travel time sensitivity kernel'
#perspective=-E175/35

ps_filename=$datafilename.ps
pdf_filename=$datafilename.pdf
jpg_filename=$datafilename.jpg

echo
echo "plotting to file: " $ps_filename
echo
verbose=-V

# ---------------------------------------------------------------------------------------------------------
# PRE-PROCESS
# ---------------------------------------------------------------------------------------------------------
# preprocess data file
if [ ! -s $datafilename ];then
  echo "can not find data file $datafilename"
  exit 1
fi

# extract data to temporary file
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

gmt gmtset PS_MEDIA letter HEADER_FONT_SIZE 14 MAP_ANNOT_OBLIQUE 0 BASEMAP_TYPE plain GMT_HISTORY false

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# without interpolation
if [ "$interpolate" = "no" ];then
  echo
  echo 'not interpolating data...'
  echo
  gmt blockmean tmpData.dat $region -I1/1 $verbose > tmpData.med.xyz
  gmt xyz2grd tmpData.med.xyz -Gimage.bin -I1/1 $region $verbose -N0.0  # -I0.5/0.5
  gmt xyz2grd tmpData.dat -Gimage.bin -I1/1 $region $verbose -N0.0      # -I0.5/0.5
  gmt grdimage image.bin $portrait $offset $Plotregion -G $projection -Bg15/g15 -T -K -V -C$colormap > $ps_filename
fi

# with interpolation
if [ "$interpolate" = "yes" ];then
  echo
  echo 'interpolating data...'
  echo
  gmt blockmean tmpData.dat $region -I1/1 $verbose > tmpData.med.xyz
  gmt surface tmpData.med.xyz  $region -I1/1 -Gimage.bin -T0.25 -C0.1 -V
  gmt grdview image.bin $perspective $zdepth -Qs -P $Plotregion $projection -K $verbose -C$colormap > $ps_filename
fi

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# store griddf  in data-file
echo "#kernel values from" $datafilename > $datafilename.xyz
echo "#lon lat val" >> $datafilename.xyz
gmt grd2xyz image.bin $verbose >> $datafilename.xyz


# image frame
gmt psscale $scalePosition $scaleAnotate -C$colormap -O $scaleoffset $portrait $verbose -K  >> $ps_filename
gmt pscoast $coastoffset $Plotregion $projection $portrait $perspective -W1 -Dc -A10000 -O -K >> $ps_filename
gmt psbasemap  $Plotregion $verbose $projection $perspective -Ba90g15/a30g15:."":WeSn -O   >> $ps_filename

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# cross-section
crosslon=45
echo '#cross-section at longitude ' $crosslon ' from' $datafilename > $datafilename.$crosslon.xyz
echo '#lon lat val' >> $datafilename.$crosslon.xyz
gawk '{if( (substr($1,1,1)!="#") && ($1!="") && (substr($1,1,1)~"[0-9]") ) \
	if( $1=='$crosslon') print($1,$2,$3); }' $datafilename.xyz >> $datafilename.$crosslon.xyz

crosslon=80
echo '#cross-section at longitude ' $crosslon ' from' $datafilename > $datafilename.$crosslon.xyz
echo '#lon lat val' >> $datafilename.$crosslon.xyz
gawk '{if( (substr($1,1,1)!="#") && ($1!="") && (substr($1,1,1)~"[0-9]") ) \
	if( $1=='$crosslon') print($1,$2,$3); }' $datafilename.xyz >> $datafilename.$crosslon.xyz

crosslon=70
echo '#cross-section at longitude ' $crosslon ' from' $datafilename > $datafilename.$crosslon.xyz
echo '#lon lat val' >> $datafilename.$crosslon.xyz
gawk '{if( (substr($1,1,1)!="#") && ($1!="") && (substr($1,1,1)~"[0-9]") ) \
	if( $1=='$crosslon') print($1,$2,$3); }' $datafilename.xyz >> $datafilename.$crosslon.xyz

crosslon=20
echo '#cross-section at longitude ' $crosslon ' from' $datafilename > $datafilename.$crosslon.xyz
echo '#lon lat val' >> $datafilename.$crosslon.xyz
gawk '{if( (substr($1,1,1)!="#") && ($1!="") && (substr($1,1,1)~"[0-9]") ) \
	if( $1=='$crosslon') print($1,$2,$3); }' $datafilename.xyz >> $datafilename.$crosslon.xyz

crosslon=10
echo '#cross-section at longitude ' $crosslon ' from' $datafilename > $datafilename.$crosslon.xyz
echo '#lon lat val' >> $datafilename.$crosslon.xyz
gawk '{if( (substr($1,1,1)!="#") && ($1!="") && (substr($1,1,1)~"[0-9]") ) \
	if( $1=='$crosslon') print($1,$2,$3); }' $datafilename.xyz >> $datafilename.$crosslon.xyz

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo
echo "converting image..."

# pdf
#convert $ps_filename $pdf_filename
#if [ -s $pdf_filename ]; then
#  rm -f $ps_filename
#fi
#open $pdf_filename

# jpg
magick -density 300 $ps_filename -quality 90 -crop 2200x1200+0+2000 $jpg_filename

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo "*******************************************************************"
echo
echo "image plotted to file           : $jpg_filename"
echo "data used is in                 : tmpData.dat"
echo "grid plotted to file            : $datafilename.xyz"
echo "cross-section at longitude 45 to: $datafilename.45.xyz"
echo
echo "*******************************************************************"



