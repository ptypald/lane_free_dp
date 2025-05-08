# Makefile for Writing Make Files Example
 
# *****************************************************
# Variables to control Makefile operation

# the compiler: gcc for C program, define as g++ for C++
CXX = clang++

# compiler flags:
#  -Wall  - this flag is used to turn on most compiler warnings
# HEADER = include/data.h
CPPFLAGS = -std=c++17 -O3 -Wall
 
main: main.o data.o dionysis_controller.o controllers.o vehicle.o idm.o dddp.o ddp.o
	$(CXX) $(CPPFLAGS) -o main main.o data.o dionysis_controller.o controllers.o vehicle.o idm.o dddp.o ddp.o

main.o: main.cpp data.h dionysis_controller.h controllers.h vehicle.h idm.h dddp.h ddp.h
	$(CXX) $(CPPFLAGS) -c main.cpp

data.o: data.h
dionysis_controller.o: dionysis_controller.h
controllers.o: controllers.h
vehicle.o: vehicle.h
idm.o: idm.h
dddp.o: dddp.h
ddp.o: ddp.h