CC = g++
CFLAGS = -pipe -O2 -std=c++14 -march=native
WARN = -Wall -Werror -Wfloat-equal -ansi -pedantic
OBJ = exchange_table.o neighbor_table.o


test.out: neighbor_table.o exchange_table.o test.o
	$(CC) $(WARN) $(CFLAGS) neighbor_table.o exchange_table.o test.o -o bin/test.out

neighbor_table.o: src/neighbor_table.cpp include/neighbor_table.h
	$(CC) $(WARN) $(CFLAGS) -c src/neighbor_table.cpp

exchange_table.o: src/exchange_table.cpp include/exchange_table.h
	$(CC) $(WARN) $(CFLAGS) -c src/exchange_table.cpp

test.o: test/test.cpp
	$(CC) $(WARN) $(CFLAGS) -c test/test.cpp

clean:
	rm -Rf *.o
