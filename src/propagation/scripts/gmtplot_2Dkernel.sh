#!/bin/bash

# ---------------------------------------------------------------------------------------------------------
# gmtplot_kernel2D.sh
#
# draws the kernel values from a given data file on a gmt map
# ---------------------------------------------------------------------------------------------------------

# arguments:
#       $1 = data-file name, e.g. ttkernel.rot.dat

# check
if [ ! -s "$1" ]; then
  echo "usage: "
  echo "    gmtplot_2Dkernel.sh kernel.dat" 
  echo
  exit
fi

# ---------------------------------------------------------------------------------------------------------
# PARAMETERS
# ---------------------------------------------------------------------------------------------------------

# interpolate grid data by gmt's spherical harmonics
interpolate=yes

# close-ups
Plotregion=-R-10/100/-60/60
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
datafilename=$1

if [ "$2" = "diff" ];then
  colormap=seis.L150.diff1.cpt
  scaleAnotate=-Ba1:"":
else
  colormap=seis.L150.4.cpt
  scaleAnotate=-Ba5:"":
fi
echo "using color map: " $colormap

title='travel time anomaly kernel - Love 150 s'
#perspective=-E175/35
ps_filename=$datafilename.ps
pdf_filename=$datafilename.pdf

echo
echo "plotting to file: " $ps_filename
echo
verbose=-V

# ---------------------------------------------------------------------------------------------------------
# PRE-PROCESS
# ---------------------------------------------------------------------------------------------------------
# preprocess data file
if [ ! -s $datafilename ];then
     echo can not find data file $datafilename
     exit
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

gmtset HEADER_FONT_SIZE 14 OBLIQUE_ANOTATION 0 BASEMAP_TYPE plain

# without interpolation
if [ "$interpolate" = "no" ];then
echo
echo 'not interpolating data...'
echo
blockmean tmpData.dat $region -I1/1 $verbose > tmpData.med.xyz
xyz2grd tmpData.med.xyz -Gimage.bin -I1/1 $region $verbose -N0.0  # -I0.5/0.5
xyz2grd tmpData.dat -Gimage.bin -I1/1 $region $verbose -N0.0      # -I0.5/0.5
grdimage image.bin $portrait $offset $Plotregion -G $projection -Bg15/g15 -T -K -V -C$colormap > $ps_filename
fi

# with interpolation
if [ "$interpolate" = "yes" ];then
echo 
echo 'interpolating data...'
echo
blockmean tmpData.dat $region -I1/1 $verbose > tmpData.med.xyz
surface tmpData.med.xyz  $region -I1/1 -Gimage.bin -T0.25 -C0.1 -V 
grdview image.bin $perspective $zdepth -Qs -P $Plotregion $projection -K $verbose -C$colormap > $ps_filename
fi
# store griddf  in data-file
echo "#kernel values from" $datafilename > $datafilename.xyz
echo "#lon lat val" >> $datafilename.xyz
grd2xyz image.bin $verbose >> $datafilename.xyz


# image frame
psscale $scalePosition $scaleAnotate -C$colormap -O $scaleoffset $portrait $verbose -K  >> $ps_filename
pscoast $coastoffset $Plotregion $projection $portrait $perspective -W1 -Dc -A10000 -O -K >> $ps_filename
psbasemap  $Plotregion $verbose $projection $perspective -Ba90g15/a30g15:."":WeSn -O   >> $ps_filename

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


echo
echo converting image to $pdf_filename

convert $ps_filename $pdf_filename
#if [ -s $pdf_filename ]; then
#  rm -f $ps_filename
#fi

echo "*******************************************************************"
echo 
echo image plotted to file: $pdf_filename
echo data used is in tmpData.dat
echo grid plotted to file: $datafilename.xyz
echo cross-section at longitude 45 to: $datafilename.45.xyz
echo
echo "*******************************************************************"

#open $pdf_filename


