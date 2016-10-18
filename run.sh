#!/bin/sh

echo 'Running Navier-stokes Solver for Tyler Green Vortex'
make -C source/

source/./result

rm source/result