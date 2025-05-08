# Copyright 2024, Gurobi Optimization, LLC

PLATFORM = macos_universal2
TWOUP    = ..
KEEPTU   = ../..
INC      = $(TWOUP)/../include/
CC       = gcc
CPP      = g++
CARGS    = -m64 -g
CPPARGS  = -std=c++98 -m64 -g
CLIB     = -L$(TWOUP)/../lib -lgurobi110
CPPLIB   = -L$(TWOUP)/../lib -lgurobi_c++ -lgurobi110
JSRC     = $(TWOUP)/java
CLASSDIR = -classpath $(TWOUP)/../lib/gurobi.jar:.
JFLAG    = -d . $(CLASSDIR)
JRUN     = -Djava.library.path=$(TWOUP)/../lib $(CLASSDIR)
GRBAPP   = dotnet
DOTNETFRAMEWORK ?= netcoreapp6.0
PROJECTFILE     ?= gurobi.csproj

all: test_c++

run: run_c++
run_c++: run_test_c++

%_c++: $(TWOUP)/c++/%_c++.cpp
	$(CPP) $(CPPARGS) -o $@ $< -I$(INC) $(CPPLIB) -lm

run_test_c++: test_c++
	./test_c++

clean:
	rm -rf *.o *_c *_c++ *.class *.log *.rlp *.lp *.bas *.ilp *.mps *.prm *.dSYM; \
	if [ -d $(GRBAPP) ]; then \
		cd $(GRBAPP); \
		find . ! \( -name "gurobi*.csproj" -o -name . \) -exec rm -rf {} +; \
	fi
