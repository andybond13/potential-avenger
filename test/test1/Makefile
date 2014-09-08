#makefile
COMPILER = g++-4
COMPILER_FLAGS = -g -O2
LIB = 
SOURCES = $(wildcard *.c) 
INC = -I .
OBJECTS = $(SOURCES:.c=.o)
EXECUTABLE = potential-avenger.exe

all: $(SOURCES) $(EXECUTABLE) 
	
$(EXECUTABLE): $(OBJECTS) 
	$(COMPILER) $(LIB) $(OBJECTS) -o $@

.c.o:
	$(COMPILER) $(COMPILER_FLAGS) $(INC) -c $< -o $@

clean:
	rm *.o *.exe *.out
