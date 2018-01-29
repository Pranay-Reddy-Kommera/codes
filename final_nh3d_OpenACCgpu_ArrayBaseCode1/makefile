.SUFFIXES:
.SUFFIXES: .f90 .o

Fxx = mpif90
FLAGS = -r8 -O3 -byteswapio -Mfree -m64 -Mnofma -Mvect=noaltcode -Minfo=accel -acc -ta=tesla:cc60
#FLAGS = -O2 -march=native
DEBUG = #-Wall -Wextra -fimplicit-none -std=f2008 -pedantic \
        #-Og -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow

OBJ = \
			basic_mod.o \
			physical_const_mod.o \
			gauss_quadrature_mod.o \
			element_mod.o \
			mpi_mod.o \
			grid_setup_mod.o  \
			testcases_nh3d_mod.o \
			prepare_edgevars_mod.o \
			dg3d_rhs_mod.o \
			euler_flux3d_mod.o \
			vertical_split_mod.o \
			timestepper_mod.o \
			gather_solution_mod.o \
			output_mod.o \
			main.o

main: $(OBJ)
	$(Fxx) $(FLAGS) $(DEBUG) -o main $(OBJ) $(LIBS)

.f90.o:
	$(Fxx) $(FLAGS) $(DEBUG) -c $*.f90

clean:
	rm -rf *.o *.mod main main.dSYM 

dclean:
	rm -rf ./data/*.dat ./out/*.out ./err/*.err
