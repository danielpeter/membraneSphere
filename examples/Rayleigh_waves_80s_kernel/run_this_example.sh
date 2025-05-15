#!/bin/bash

# setup bin/ and data/ folder in this example directory
if [ ! -e bin ]; then ln -s ../../bin/; fi
if [ ! -e data ]; then ln -s ../../data/; fi

# clean
mkdir -p OUTPUT; rm -rf OUTPUT/*

# backup
cp -vp Parameter_Input Parameter_Input.org

##
## combined forward/adjoint  - simulation
##
echo
echo "***********************************************"
echo "running adjoint method..."
echo "***********************************************"
echo

# run simulation
mpirun -np 2 ./bin/adjointMethod | tee OUTPUT/output_adjointMethod.txt

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo
echo "done"
echo

