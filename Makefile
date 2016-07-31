#Makefile for Liquid N2 fitting
#Simon Michalke, (c) 2016

CC=g++
#CFLAGS=-pthread -std=c++11 -fPIC -c -O3 #optimized
CFLAGS=-pthread -std=c++11 -fPIC -c -Wall -g #debugging
LDFLAGS=-pthread -std=c++11
INCLDIR=

all: main

clean:
	@touch temp/dummy
	rm temp/*

main: temp/main.o temp/vector.o temp/simulation.o
	$(CC) -g -o liquid_n2_fitting \
	temp/main.o temp/vector.o temp/simulation.o \
	$(LDFLAGS)

#all source files:
temp/main.o: src/main.cpp src/define.h src/simulation.h temp
	$(CC) $(CFLAGS) $(INCLDIR) src/main.cpp   -o temp/main.o
temp/vector.o: src/vector.cpp src/vector.h temp
	$(CC) $(CFLAGS) $(INCLDIR) src/vector.cpp -o temp/vector.o
temp/simulation.o: src/simulation.cpp src/simulation.h temp
	$(CC) $(CFLAGS) $(INCLDIR) src/simulation.cpp -o temp/simulation.o

#header files:
src/define.h:
src/vector.h:     src/define.h
src/simulation.h: src/define.h src/vector.h

#@touch temp/dummy
temp:
	@mkdir -p temp
