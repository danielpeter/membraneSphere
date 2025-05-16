#!/bin/bash

level=$1
if [ "$level" == "" ]; then echo "usage: ./plot_grid.sh level[e.g., 6]"; exit 1; fi

# setup script file
if [ ! -e ./plot_3Dgrid.py ]; then
  ln -s ../../scripts/plot_3Dgrid.py
fi

echo
echo "***********************************************"
echo "plotting grid cellFractionAverages ..."
echo "***********************************************"
echo

./plot_3Dgrid.py OUTPUT/Dvvert${level}.dat OUTPUT/Dvface${level}.dat OUTPUT/cellFractionAverage.${level}.dat

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

mv -v out.jpg OUTPUT/cellFractionAverage.${level}.dat.jpg

echo
echo "done"
echo

