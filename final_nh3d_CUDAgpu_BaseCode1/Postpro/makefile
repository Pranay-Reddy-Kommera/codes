.SUFFIXES:
.SUFFIXES: .f90 .o

Fxx = mpif90
FLAGS = -O2
#FLAGS = -O2 -march=native
DEBUG = #-Wall -Wextra -fimplicit-none -std=f2008 -pedantic \
        #-Og -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow

OBJ = \
			basic_mod.o \
			gauss_quadrature_mod.o \
			interpol_mod.o \
			grids_mod.o \
			post_process.o

main: $(OBJ)
	$(Fxx) $(FLAGS) $(DEBUG) -o postp $(OBJ) $(LIBS)

.f90.o:
	$(Fxx) $(FLAGS) $(DEBUG) -c $*.f90

clean:
	rm -rf *.o *.mod postp 

dclean:
	rm -rf ./data/*.dat
