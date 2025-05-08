# Makefile for Writing Make Files Example
 
# *****************************************************
# Variables to control Makefile operation

# the compiler: gcc for C program, define as g++ for C++
CC = clang++

# compiler flags:
#  -Wall  - this flag is used to turn on most compiler warnings
CFLAGS = -Wall -std=c++11 -O3
 
TARGET = main

all: $(TARGET)

$(TARGET): $(TARGET).cpp
	$(CC) $(CFLAGS) -o $(TARGET) $(TARGET).cpp

clean:
	$(RM) $(TARGET)