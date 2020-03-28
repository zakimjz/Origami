### GLOBALS
#CC	   	= g++ -DPRINT
CC	   	= g++
BOOSTLIB        =/usr/local/lib
CFLAGS	   	= -g -O3 -Wall -L$(BOOSTLIB)
# CFLAGS	   	= -g -O3 -Wall -include /usr/local/include/mpatrol.h -lmpatrol -lbfd
INCLUDE-PATH 	= -I. -I../src/common -I../src/graph -I../src/StringTokenizer
                  
INCLUDES-COMMON = ../src/common/*.h
INCLUDES-GRAPH = ../src/graph/*.h graph_test.cpp
INCLUDES-TOKEN = ../src/StringTokenizer/*.h
OBJ            = ../src/StringTokenizer/StringTokenizer.o
### TARGETS
MEMORY-BASED  = graph_test

all: 
	cd ../src/StringTokenizer; 	$(MAKE);
	$(MAKE) $(MEMORY-BASED)

memory-based: $(MEMORY-BASED)

default: all

clean:
	cd ../src/StringTokenizer; 	$(MAKE) clean;
	rm -f $(MEMORY-BASED) *.o 

### RULES
.SUFFIXES: .cpp

.cpp: 
	$(CC) $(CFLAGS) $(INCLUDE-PATH) $(OBJ) $< -o $@

### DEPENDENCIES
# target: headers
graph_test:  $(INCLUDES-GRAPH) $(INCLUDES-COMMON)
