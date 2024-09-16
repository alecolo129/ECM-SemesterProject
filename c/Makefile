CC = gcc
all: test

clean: 
	rm -f ecm test *.o

test: test.c ecm.h curves.h
	$(CC) -O3 -o test test.c -lgmp 

.PHONY: all clean

