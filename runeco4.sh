#
# Runscript for ecosys_test:
#
# ------------------ COMPILE & RUN ----------------------------
#
 gfortran -fbounds-check  mod_geochem.F90 mod_param.F90 mod_heat.F90 mod_foodweb3.F90 mod_coral3.F90 mod_algae.F90 mod_seagrass.F90 mod_sedecosys.F90 mod_reef_ecosys3.F90 mod_input.F90 ecosys_test5.F90 -O2 -o ecosys_test5.exe
#
# ------------------ OUTPUT FILES -----------------------------
#
 ./ecosys_test5.exe
#
