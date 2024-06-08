CC = g++
CFLAGS = -O3 --std=c++11
all: es07.2.x
es07.2.x: es07.2.o
	$(CC) es07.2.o -o es07.2.x
es07.2.o: es07.2.cpp
	$(CC) -c es07.2.cpp -o es07.2.o $(CFLAGS)
clean :
	rm *.o es07.2.x