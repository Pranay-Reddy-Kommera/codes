# Makefile

#CC = pgc++
CC = icc
LIB = -L/home/sureshm/pranay/softwares/SuiteSparse/lib -lcholmod -lamd -lcolamd -lsuitesparseconfig

CFLAGS = -std=c++11 -I/home/sureshm/pranay/softwares/SuiteSparse/include
DEP = BA_includeFiles.h BA_data.h
OBJ = BA_main.o BA_dataRW.o BA_functions.o BA_subFunctions.o BA_Math.o

%.o: %.cpp $(DEP)
	$(CC) -c -o $@ $< $(CFLAGS) $(LIB)

ba_schur: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIB)

clean:
	rm -rf *.o ba_schur
