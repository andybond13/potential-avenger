#makefile
COMPILER = g++-5
COMPILER_FLAGS = -g -O3 -Wall -std=c++11 -fopenmp -fcilkplus  -DCILKPLUS_GCC=1 #-ftree-parallelize-loops=64
LIB = -L/sw/opt/boost-1_55/lib/ -lboost_random -lcilkrts #-fopenacc 
SOURCES = $(wildcard *.c) 
INC = -I . -I/sw/opt/boost-1_55/include/  #-fopenacc
OBJECTS = $(SOURCES:.c=.o)
EXECUTABLE = potential-avenger.exe

all: $(SOURCES) $(EXECUTABLE) dSYM
	
$(EXECUTABLE): $(OBJECTS) 
	$(COMPILER) $(LIB) $(OBJECTS) -o $@

.c.o:
	$(COMPILER) $(COMPILER_FLAGS) $(INC) -c $< -o $@

dSYM:
	dsymutil potential-avenger.exe -o potential-avenger.exe.dSYM

runtest: clean all
	cd test && python makeAndRunTests.py

clean:
	rm *.o *.exe; rm -rf *dSYM
