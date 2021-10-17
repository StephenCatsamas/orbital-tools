SRC = main.cpp util.cpp integrator.cpp
INC = main.h util.h integrator.h 
CFLAGS = -Wall 
OBJ = orbits.exe

all: $(SRC) $(INC)
	g++ $(CFLAGS) $(SRC) -o $(OBJ)
	
run: all
	orbits.exe && vis.py
	
profile:
	g++ -g -pg -no-pie -fno-builtin main.cpp -o out.exe
	

debug: $(SRC) $(INC)
	g++ $(CFLAGS) $(SRC) -g
	gdb a.exe
