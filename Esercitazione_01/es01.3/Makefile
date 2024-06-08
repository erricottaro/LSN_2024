CC = g++
CFLAGS = -O3 --std=c++11
HEADERS = random.h

all: es01.3.x test.x

test.x: test.o random.o
	$(CC) random.o test.o -o test.x $(CFLAGS)
es01.3.x: es01.3.o random.o
	$(CC) random.o es01.3.o -o es01.3.x $(CFLAGS)
%.o: %.cpp $(HEADERS)
	$(CC) -c $< -o $@ $(CFLAGS)

clean :
	rm *.o *.x
