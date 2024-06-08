CC = g++
CFLAGS = -O3 --std=c++11
HEADERS = random.h

es01.2.x: es01.2.o random.o
	$(CC) random.o es01.2.o -o es01.2.x $(CFLAGS)
%.o: %.cpp $(HEADERS)
	$(CC) -c $< -o $@ $(CFLAGS)

clean :
	rm *.o *.x
