SRC = main.cpp util.cpp integrator.cpp
INC = main.h util.h integrator.h 
CFLAGS = -Wall -O2 -ftree-vectorize
OBJ = orbits.exe
CC = g++

all: $(SRC) $(INC)
	$(CC) $(CFLAGS) $(SRC) -o $(OBJ)
	
run: all
	orbits.exe && vis.py
	
profile: $(SRC) $(INC)
	$(CC) -g -gdwarf-2 -shared-libgcc -D_DEBUG $(CFLAGS) $(SRC) -o $(OBJ)
#-fno-pie -fno-builtin
	
debug: $(SRC) $(INC)
	$(CC) $(CFLAGS) $(SRC) -g
	gdb a.exe