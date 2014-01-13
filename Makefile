#makefile
COMPILER = g++
COMPILER_FLAGS = -g -O3 -Wall
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
	rm *.o *.exe
