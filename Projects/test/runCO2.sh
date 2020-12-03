#
# Runscript for ecosys_test:
#
# ------------------ COMPILE & RUN ----------------------------
#
#ifort pom08.f -o a.out /usr/local/lib/libnetcdf.a
 gfortran mod_geochem.F90 CO2system_test.f90 -o CO2system_test.exe
#
# ------------------ OUTPUT FILES -----------------------------
#
#a.out > pom08.out        # printout file
 ./CO2system_test.exe
#
#