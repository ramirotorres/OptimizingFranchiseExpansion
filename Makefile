# Copyright 2020, Gurobi Optimization, LLC

PLATFORM = mac64
PATH    = /Library/gurobi903/mac64
INC      = $(PATH)/include/
CPP      = g++
CARGS    = -m64 -g
CPPLIB   = -L$(PATH)/lib -lgurobi_c++ -lgurobi90

all: Modelo Modelo_demo Modelo_Heuristica

Modelo:Modelo.cpp 
	$(CPP) $(CARGS)  -o $@ $< -I$(INC) $(CPPLIB) -lm

Modelo_demo:Modelo_demo.cpp 
	$(CPP) $(CARGS)  -o $@ $< -I$(INC) $(CPPLIB) -lm

Modelo_Heuristica:Modelo_Heuristica.cpp 
	$(CPP) $(CARGS)  -o $@ $< -I$(INC) $(CPPLIB) -lm

clean:
	rm -rf *.o;
