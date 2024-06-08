CC = g++
CFLAGS = -O3 --std=c++11
HEADERS = random.h

all: es02.2.x es02.2.spherical.x

es02.2.spherical.x: es02.2.spherical.o random.o 
	$(CC) random.o es02.2.spherical.o -o es02.2.spherical.x $(CFLAGS)
es02.2.x: es02.2.o random.o
	$(CC) random.o es02.2.o -o es02.2.x $(CFLAGS)
%.o: %.cpp $(HEADERS)
	$(CC) -c $< -o $@ $(CFLAGS)

clean :
	rm *.o *.x
