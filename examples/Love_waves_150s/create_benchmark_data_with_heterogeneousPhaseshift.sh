#!/bin/bash
#
# reads in (real) phase shift data (e.g., wei_sum.02.L0150.1.txt) and uses the same source/receiver data setups
# to create a new data file for a given heterogeneouse phase map, which can be used for a benchmark inversion
#
# setup bin/ and data/ folder in this example directory
if [ ! -e bin ]; then ln -s ../../bin/; fi
if [ ! -e data ]; then ln -s ../../data/; fi

# clean
mkdir -p OUTPUT; rm -rf OUTPUT/*

# backup
cp -vp Parameter_Input Parameter_Input.org

##
## heterogeneousePhaseshift
##
echo
echo "***********************************************"
echo "running heterogeneousePhaseshift ..."
echo "***********************************************"
echo

# full dataset or only a small subset test
TEST_ONLY=1

# only for testing
if [ "$TEST_ONLY" == "1" ]; then
  # create a small subset only for testing
  head -n 25 data/phasedata/wei_sum.02.L0150.1.txt > wei_sum.short.txt
  # updates data file
  sed "s:^INV_DATA .*:INV_DATA                         = 'wei_sum.short.txt':" Parameter_Input.org > Parameter_Input
fi

# run binary
mpirun -np 8 ./bin/heterogeneousPhaseshift | tee OUTPUT/output_heterogeneousPhaseshift.txt

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo
echo "done"
echo

