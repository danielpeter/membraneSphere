#!/bin/bash

# setup bin/ and data/ folder in this example directory
if [ ! -e bin ]; then ln -s ../../bin/; fi
if [ ! -e data ]; then ln -s ../../data/; fi

# clean
mkdir -p OUTPUT; rm -rf OUTPUT/*

# backup
cp -vp Parameter_Input Parameter_Input.org

##
## homogeneous phase map - simulation
##
echo
echo "***********************************************"
echo "running propagation w/ homogeneous phase map..."
echo "***********************************************"
echo
sed "s:^HETEROGENEOUS                    =:HETEROGENEOUS                    = .false.:" Parameter_Input.org > Parameter_Input

# run simulation
./bin/propagation | tee OUTPUT/output_propagation.hom.txt

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# rename seismogram output
mv -v OUTPUT/seismo.L150.withoutDelta.dat OUTPUT/seismo.hom.dat


##
## heterogeneous phase map - simulation
##
echo
echo "***********************************************"
echo "running propagation w/ heterogeneous phase map..."
echo "***********************************************"
echo
sed "s:^HETEROGENEOUS                    =:HETEROGENEOUS                    = .true.:" Parameter_Input.org > Parameter_Input

# run simulation
./bin/propagation | tee OUTPUT/output_propagation.het.txt

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# rename seismogram output
mv -v OUTPUT/seismo.L150.withoutDelta.dat OUTPUT/seismo.het.dat

##
## timelag
##
echo
echo "***********************************************"
echo "running timelag..."
echo "***********************************************"
echo
./bin/timelag | tee OUTPUT/output_timlag.txt

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo
echo "done"
echo

