CC = g++
override CFLAGS += -Wall -lgmp -O3 -I /Users/depo/Documents/Uni/4_anno/SemesterProject/Re-discoveringECM/c++/gmp-6.2.1 -L /Users/depo/Documents/Uni/4_anno/SemesterProject/Re-discoveringECM/c++/gmp-6.2.1 -std=c++17 

all: ecm old

clean: 
	rm -f ecm old *o

.PHONY: all clean

ecm: parameters.h ecm.cpp
	$(CC) $(CFLAGS) -o ecm ecm.cpp

old: parameters.h old.cpp
	$(CC) $(CFLAGS) -o old old.cpp

