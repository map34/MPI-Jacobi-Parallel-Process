CFLAGS=-g -Werror=vla
AFTER_FLAG=-lm

all: main.exe

main.exe: main.o MyFunctions.o
	mpicc $(CFLAGS) -o main.exe main.o MyFunctions.o $(AFTER_FLAG)
	make clean_o

main.o: main.c MyHeader.h
	mpicc $(CFLAGS) -c main.c $(AFTER_FLAG)

MyFunctions.o: MyFunctions.c MyHeader.h
	mpicc $(CFLAGS) -c MyFunctions.c $(AFTER_FLAG)

clean:
	rm *.o main

clean_o:
	rm *.o
