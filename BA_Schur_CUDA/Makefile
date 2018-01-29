# Makefile

#CC = pgc++
CC = nvcc
LIB = -L/usr/local/cuda-9.0/lib64 -lcudart -lcusparse

CFLAGS = -std=c++11 -arch=sm_35 -I/usr/local/cuda-9.0/include
DEP = BA_includeFiles.h BA_data.h
OBJ = BA_main.o BA_dataRW.o BA_functions.o BA_subFunctions.o BA_Math.o BA_cudaMath.o

%.o: %.cpp $(DEP)
	$(CC) -c -o $@ $< $(CFLAGS) $(LIB)

%.o: %.cu $(DEP)
	$(CC) -c -o $@ $< $(CFLAGS) $(LIB)

ba_schur: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIB)

clean:
	rm -rf *.o ba_schur
