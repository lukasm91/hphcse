# Makefile for n-particle problem
CC = g++
CFLAGS = -Wall -pedantic -ansi -c
LFLAGS = -Wall
nBody: main.o domain.o particle.o
	$(CC) $(LFLAGS) -o $@ $^
main.o: main.cpp domain.h
	$(CC) $(CFLAGS) $<
domain.o: domain.cpp domain.h particle.h
	$(CC) $(CFLAGS) $<
particle.o: particle.cpp particle.h
	$(CC) $(CFLAGS) $<

