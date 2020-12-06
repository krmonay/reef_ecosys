#!/bin/sh
#
SRC_DIR=../../src
INCLUDE="-I${PWD}"
gfortran -fbounds-check -ffree-form ${SRC_DIR}/mod_geochem.F ${SRC_DIR}/mod_param.F ${SRC_DIR}/mod_foodweb.F ${SRC_DIR}/mod_coral.F ${SRC_DIR}/mod_macroalgae.F ${SRC_DIR}/mod_seagrass.F ${SRC_DIR}/mod_sedecosys.F ${SRC_DIR}/mod_reef_ecosys.F ${SRC_DIR}/mod_input.F ${SRC_DIR}/mod_output.F ${SRC_DIR}/ecosys_test.F ${INCLUDE} -O2 -o ecosys_test.exe
rm *.mod
#
mkdir -p output
#
./ecosys_test.exe
#
