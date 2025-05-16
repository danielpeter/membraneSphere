#!/bin/bash

level=$1
if [ "$level" == "" ]; then echo "usage: ./plot_grid.sh level[e.g., 6]"; exit 1; fi

# setup script file
if [ ! -e ./plot_3Dgrid.py ]; then
  ln -s ../../scripts/plot_3Dgrid.py
fi

echo
echo "***********************************************"
echo "plotting numerical approximation error ..."
echo "***********************************************"
echo

./plot_3Dgrid.py data/griddata/Dvvert${level}.dat data/griddata/Dvface${level}.dat OUTPUT/numericalError${level}.L6M1.dat

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

mv -v out.jpg OUTPUT/numericalError${level}.L6M1.dat.jpg

echo
echo "done"
echo

