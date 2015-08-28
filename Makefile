#makefile
COMPILER_FLAGS = -g -Wall
LIB = -L/sw/opt/boost-1_55/lib/ -lboost_random 
COMPILER = g++-5
SOURCES = $(wildcard *.c) 
INC = -I . -I/sw/opt/boost-1_55/include/
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
