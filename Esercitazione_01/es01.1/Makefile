CC = g++
CFLAGS = -O3 --std=c++11
HEADERS = random.h

all: es01.1.x es01.1chisquared.x

es01.1chisquared.x: es01.1chisquared.o random.o
	$(CC) random.o es01.1chisquared.o -o es01.1chisquared.x $(CFLAGS)

es01.1.x: es01.1.o random.o
	$(CC) random.o es01.1.o -o es01.1.x $(CFLAGS)

%.o: %.cpp $(HEADERS)
	$(CC) -c $< -o $@ $(CFLAGS)

clean :
	rm *.o *.x
