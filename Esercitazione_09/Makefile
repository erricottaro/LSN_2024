CC = g++
CFLAGS = -O3 --std=c++11
HEADERS = random.h city.h individual.h population.h
all: es09.1
es09.1: es09.1.o individual.o city.o random.o population.o
	$(CC) population.o individual.o city.o random.o es09.1.o -o es09.1
%.o: %.cpp $(HEADERS)
	$(CC) -c $< -o $@ $(CFLAGS)
clean :
	rm *.o *.x
remove:
	rm OUTPUT/*.*
