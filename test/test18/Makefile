#makefile
COMPILER = g++-5
COMPILER_FLAGS = -O2 -Wall
LIB = -L/sw/opt/boost-1_55/lib/ -lboost_random -fopenmp -lcilkrts 
SOURCES = $(wildcard *.c) 
INC = -I . -I/sw/opt/boost-1_55/include/ -fopenmp -fcilkplus
OBJECTS = $(SOURCES:.c=.o)
EXECUTABLE = potential-avenger.exe

all: $(SOURCES) $(EXECUTABLE) 
	
$(EXECUTABLE): $(OBJECTS) 
	$(COMPILER) $(LIB) $(OBJECTS) -o $@

.c.o:
	$(COMPILER) $(COMPILER_FLAGS) $(INC) -c $< -o $@

clean:
	rm *.o *.exe *.out
