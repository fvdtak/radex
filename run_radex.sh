#!/bin/bash
#export OMP_STACKSIZE=512M
#ulimit -s unlimited
cd src
make clean
make
cd ../
#executable and input_params
./radexx < test_hco_var1.inp
