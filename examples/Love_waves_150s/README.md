# membrane waves - Love waves at 150 s

This example simulates Love waves at 150 s period.
It performs two simulations, one for a homogeneous phase-velocity map,
and a second one with a heterogeneous phase-velocity map.

The simulation input parameters are in file `Parameter_Input`, the wave propagation executable is `bin/propagation`.
After the two simulations, the cross-correlation time lag is computed between the two resulting traces
by the executable `bin/timelag`, which uses the input parameters from file `Timelag_Input`.

To run this example, type:
```
$ ./run_this_example.sh
```


