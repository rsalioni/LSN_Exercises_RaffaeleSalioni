CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o Integral.o Funzioni.o
	$(CC) random.o main.o Integral.o Funzioni.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
Integral.o : Integral.cpp Integral.h
	$(CC) -c Integral.cpp -o Integral.o $(CFLAGS)
Funzioni.o : Funzioni.cpp Funzioni.h
	$(CC) -c Funzioni.cpp -o Funzioni.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out
