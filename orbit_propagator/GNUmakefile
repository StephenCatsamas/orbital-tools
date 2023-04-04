SRC = main.cpp util.cpp integrator.cpp interface.cpp
INC = main.h util.h integrator.h 
CFLAGS = -Wall -Wextra -Wno-unused-result -Wno-format-extra-args \
         -O2 -ftree-vectorize
# CFLAGS = -Wall
OBJ = orbits
CC = g++

all: $(SRC) $(INC)
	$(CC) $(CFLAGS) $(SRC) -o $(OBJ)
	
run: all
	./$(OBJ) && py ./vis.py
	
profile: $(SRC) $(INC)
	$(CC) -g -gdwarf-2 -shared-libgcc -D_DEBUG $(CFLAGS) $(SRC) -o $(OBJ)

	
debug: $(SRC) $(INC)
	$(CC) $(CFLAGS) $(SRC) -g
	gdb a.exe
