CC = g++
CFLAGS = -O3 --std=c++11
HEADERS = random.h

all: es02.1.x

es02.1.x: es02.1.o random.o
	$(CC) random.o es02.1.o -o es02.1.x $(CFLAGS)
%.o: %.cpp $(HEADERS)
	$(CC) -c $< -o $@ $(CFLAGS)

clean :
	rm *.o *.x
