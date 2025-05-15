#!/bin/bash

level=$1
if [ "$level" == "" ]; then echo "usage: ./run_grid_construction_level.sh level[e.g., 6]"; exit 1; fi

# setup bin/ and data/ folder in this example directory
if [ ! -e bin ]; then ln -s ../../bin/; fi
if [ ! -e data ]; then ln -s ../../data/; fi

# check binaries
if [ ! -e bin/gridConstruction ]; then
  echo "Please compile first all grid binaries in the root directory by 'make grid', before running this example..."; echo
  exit 1
  #cd ../../
  #make grid
  #cd $currentdir
fi

# clean
mkdir -p OUTPUT; rm -rf OUTPUT/*

echo `date`
echo "***********************************************"
echo "running grid construction up to level $level ..."
echo "***********************************************"
echo

./bin/gridConstruction <<EOF
$level
EOF

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# transfer?
echo
echo "(optional) to transfer grid data:"
echo "  > rsync -avz ./OUTPUT/D*.dat ./data/griddata/"
echo

echo
echo `date`
echo "done"
echo

