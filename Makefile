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

main: temp/main.o
	$(CC) -g -o liquid_n2_fitting \
	temp/main.o \
	$(LDFLAGS)

#all source files:
temp/main.o: src/main.cpp temp
	$(CC) $(CFLAGS) $(INCLDIR) src/main.cpp   -o temp/main.o

#header files:

#@touch temp/dummy
temp:
	@mkdir -p temp
