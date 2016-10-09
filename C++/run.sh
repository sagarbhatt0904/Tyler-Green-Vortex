#!/bin/sh

echo 'Running Navier-stokes Solver for Tyler Green Vortex'
make -C source/

source/./tyler_green

rm source/tyler_green