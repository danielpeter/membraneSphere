#!/bin/bash

level=$1
if [ "$level" == "" ]; then echo "usage: ./run_grid_Laplacian_numericalError.sh level[e.g., 6]"; exit 1; fi

# setup bin/ and data/ folder in this example directory
if [ ! -e bin ]; then ln -s ../../bin/; fi
if [ ! -e data ]; then ln -s ../../data/; fi

# check binaries
if [ ! -e bin/sphericalLaplacian ]; then
  echo "Please compile first all grid binaries in the root directory by 'make auxiliaries', before running this example..."; echo
  exit 1
  #cd ../../
  #make auxiliaries
  #cd $currentdir
fi

# output directory
mkdir -p OUTPUT;

echo
echo "***********************************************"
echo "running Laplacian accuracy ..."
echo "***********************************************"
echo

./bin/sphericalLaplacian <<EOF
$level
EOF

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo
echo "done"
echo

