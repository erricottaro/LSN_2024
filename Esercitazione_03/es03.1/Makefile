CC = g++
CFLAGS = -O3 --std=c++11
HEADERS = random.h

all: es03.1.x

es03.1.x: es03.1.o random.o
	$(CC) random.o es03.1.o -o es03.1.x $(CFLAGS)
%.o: %.cpp $(HEADERS)
	$(CC) -c $< -o $@ $(CFLAGS)

clean :
	rm *.o *.x
