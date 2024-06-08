CC = mpicxx
CFLAGS = -O3 --std=c++11
HEADERS = random.h city.h individual.h population.h
all: es10.2
es10.2: es10.2.o individual.o city.o random.o population.o
	$(CC) population.o individual.o city.o random.o es10.2.o -o es10.2
%.o: %.cpp $(HEADERS)
	$(CC) -c $< -o $@ $(CFLAGS)
clean :
	rm *.o *.x
remove:
	rm OUTPUT/*.*
