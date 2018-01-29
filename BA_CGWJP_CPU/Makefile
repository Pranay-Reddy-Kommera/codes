# Makefile

#CC = pgc++
CC = icc
#LIB = -L/home/sureshm/pranay/softwares/SuiteSparse/lib -lcholmod -lamd -lcolamd -lsuitesparseconfig

CFLAGS = -std=c++11
DEP = BA_include.h BA_declaration.h
OBJ = BA_main.o BA_io.o BA_math.o BA_functions.o BA_subFunctions.o

%.o: %.cpp $(DEP)
	$(CC) -c -o $@ $< $(CFLAGS) $(LIB)

ba_cgwop: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIB)

clean:
	rm -rf *.o ba_cgwop
