#!/bin/bash

level=$1
if [ "$level" == "" ]; then echo "usage: ./plot_grid.sh level[e.g., 6]"; exit 1; fi

# setup script file
if [ ! -e ./plot_3Dgrid.py ]; then
  ln -s ../../scripts/plot_3Dgrid.py
fi

echo
echo "***********************************************"
echo "plotting spherical harmonics ..."
echo "***********************************************"
echo

./plot_3Dgrid.py data/griddata/Dtvert${level}.dat data/griddata/Dtface${level}.dat OUTPUT/harmonics${level}.L6M1.dat stretch=0.5

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

mv -v out.jpg OUTPUT/harmonics${level}.L6M1.dat.jpg

echo
echo "done"
echo

