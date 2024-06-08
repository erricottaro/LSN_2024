CC = g++
CFLAGS = -O3 --std=c++11
all: es08.1.x es08.2.x
es08.2.x: annealing.o es08.2.o system.o random.o
	$(CC) annealing.o system.o random.o es08.2.o -o es08.2.x
es08.2.o: es08.2.cpp annealing.h
	$(CC) -c es08.2.cpp -o es08.2.o $(CFLAGS) 
es08.1.x : system.o es08.1.o random.o
	$(CC) system.o random.o es08.1.o -o es08.1.x
es08.1.o : es08.1.cpp system.h
	$(CC) -c es08.1.cpp -o es08.1.o $(CFLAGS)
annealing.o: annealing.cpp system.h random.h
	$(CC) -c annealing.cpp -o annealing.o
system.o : system.cpp system.h random.h
	$(CC) -c system.cpp -o system.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o *.x
remove:
	rm OUTPUT/*.*
