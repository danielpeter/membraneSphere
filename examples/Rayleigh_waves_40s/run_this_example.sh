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

latitudes=(90 60 30 0 -30 -60 -90)
for lat in "${latitudes[@]}"; do
  # Format the integer into a float with 1 decimal place
  f_lat=$(printf "%.1f" "$lat")
  echo
  echo "receiver latitude: $f_lat"
  echo 
  # update Parameter files
  sed "s/RECEIVER .*/RECEIVER                         = ${f_lat} 0.0/" Parameter_Input.org > Parameter_Input

  # run simulation
  mpirun -np 1 ./bin/propagation | tee OUTPUT/output_propagation.txt

  # checks exit code
  if [[ $? -ne 0 ]]; then exit 1; fi

  # rename outputs
  mv -v OUTPUT/output_propagation.txt OUTPUT/output_propagation.lat${lat}.txt
  mv -v OUTPUT/seismo.R40.withoutDelta.dat OUTPUT/seismo.R40.withoutDelta.lat${lat}.dat
done

echo
echo "done"
echo

