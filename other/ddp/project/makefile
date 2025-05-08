# Makefile for Writing Make Files Example
 
# *****************************************************
# Variables to control Makefile operation

# the compiler: gcc for C program, define as g++ for C++
CC = clang++

# compiler flags:
#  -Wall  - this flag is used to turn on most compiler warnings
CFLAGS = -Wall -std=c++17 -O3
 
main: main.o data.o dionysis_controller.o controllers.o vehicle.o
	$(CC) $(CFLAGS) -o main main.o data.o dionysis_controller.o controllers.o vehicle.o

main.o: main.cpp data.h dionysis_controller.h controllers.h vehicle.h
	$(CC) $(CFLAGS) -c main.cpp

data.o: data.h
dionysis_controller.o: dionysis_controller.h
controllers.o: controllers.h
vehicle.o: vehicle.h
