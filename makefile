# Name of executable
PROGRAM=hopping

# List of all object files in order of compilation
OBJ=types.o maths.o system.o fmo.o sbm.o pyr.o dma.o ful.o so2.o tul.o rhe.o dynamics.o

# External libraries
LIBS=-llapack

# Choose between "debug" and "fast" compilation
default: debug

# Debug checks for plenty of bugs
debug: FC=gfortran -g -fPIC -Wall -Wextra -Wconversion -pedantic --check=all -Wuninitialized -fbacktrace -finit-real=nan -fbounds-check #-ffpe-trap=invalid,zero,overflow
debug: main

# Fast enables parallelization
fast: FC=gfortran -O2 -fopenmp -fPIC
fast: main

# Compile main script
main: $(OBJ) main.f90
	$(FC) -o $(PROGRAM).x $^ $(LIBS)

# Compiles modules whose source files have changed since last compilation
%.o: %.f90
	$(FC) -c $^

%.o: %.f
	$(FC) -c $^

# Run "make clean" to clean up. Do this when switching between debug and fast compilation
clean:
	rm -f -r *.o *.x *.mod *.x.dSYM
