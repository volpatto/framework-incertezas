# Clearing
binDir="bin/" ; comando="cd $binDir"; eval $comando
comando="rm *.o"; eval $comando
comando="rm *.mod"; eval $comando

# Compiling modules
comando="gfortran -c ../src/random.f90"; eval $comando
comando="gfortran -c ../src/utilities.f90"; eval $comando
comando="gfortran -c ../src/tintegration.f90"; eval $comando
comando="gfortran -c ../src/uncertainty.f90"; eval $comando
# Linking and compiling driver
#comando="gfortran ../src/main.f90 asa183.o ranlib.o utilities.o uncertainty.o -o ../main.exe"; eval $comando
#comando="gfortran ../src/main.f90 asa183.o nrtype.o nrutil.o ran_state.o ranrecipies.o utilities.o tintegration.o uncertainty.o -o ../main.exe"; eval $comando
comando="gfortran ../src/main.f90 random.o utilities.o tintegration.o uncertainty.o -o ../main.exe"; eval $comando
#comando="gfortran ../src/main.f90 utilities.o uncertainty.o -o ../main.exe"; eval $comando
