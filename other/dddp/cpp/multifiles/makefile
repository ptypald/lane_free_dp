# Makefile for Writing Make Files Example
 
# *****************************************************
# Variables to control Makefile operation

# the compiler: gcc for C program, define as g++ for C++
CC = clang++

# compiler flags:
#  -Wall  - this flag is used to turn on most compiler warnings
CFLAGS = -Wall -std=c++11 -O3
 
main: main.o data.o vehicle.o controller.o dddp.o
	$(CC) $(CFLAGS) -o main main.o data.o vehicle.o controller.o dddp.o

main.o: main.cpp data.h infrastructure.h vehicle.h controller.h dddp.h
	$(CC) $(CFLAGS) -c main.cpp

data.o: data.h
vehicle.o: vehicle.h
controller.o: controller.h
dddp.o: dddp.h

# ~~~~~~~~~~~~~~~~~~