CC = g++
CFLAGS = -pipe -O2 -std=c++14 -march=native
WARN = -Wall -Werror -Wfloat-equal -ansi -pedantic
OBJ = exchange_table.o neighbor_table.o

test_neigh.out: neighbor_table.o test_neigh.o
	$(CC) $(WARN) $(CFLAGS) neighbor_table.o test_neigh.o -o test_neigh.out

neighbor_table.o: src/neighbor_table.cpp include/neighbor_table.h
	$(CC) $(WARN) $(CFLAGS) -c src/neighbor_table.cpp

test_neigh.o: test/test_neigh.cpp
	$(CC) $(WARN) $(CFLAGS) -c test/test_neigh.cpp

