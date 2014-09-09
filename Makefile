#makefile
COMPILER = g++-4
COMPILER_FLAGS = -g -Wall
LIB = 
SOURCES = $(wildcard *.c) 
INC = -I .
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
