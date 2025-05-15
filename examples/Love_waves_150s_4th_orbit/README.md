# membrane waves - Love waves at 150 s with 4th-orbit comparison

This example simulates Love waves at 150 s period up to the 4th-orbit.
This is similar to the example in folder `Love_waves_150s/`,
but it extends the time span of the simulation
and runs in parallel on 4 MPI processes.

It performs two simulations, one for a homogeneous phase-velocity map,
and a second one with a heterogeneous phase-velocity map.
After the two simulations, the cross-correlation time lag is computed between the two resulting traces,
focussing only on the 4th-orbit signal differences.

To run this example, type:
```
$ ./run_this_example.sh
```


