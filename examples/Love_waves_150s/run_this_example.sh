#!/bin/bash

# setup bin/ and data/ folder in this example directory
if [ ! -e bin ]; then
  ln -s ../../bin/
fi
if [ ! -e data ]; then
  ln -s ../../data/
fi

mkdir -p OUTPUT

# run simulation
./bin/propagation | tee OUTPUT/output_propagation.txt


echo
echo "done"
echo

