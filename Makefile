SHELL = /bin/sh
CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -funroll-loops -std=c++0x -O3 -DNDEBUG -march=native -pthread
#CFLAGS = -Wall -Wextra -pedantic -Wshadow -std=c++0x -g2 -pthread

all: decode_lz77
decode_lz77:
	$(CC) $(CFLAGS) -o decode_lz77 src/main.cpp src/utils.cpp
clean:
	/bin/rm -f *.o
nuclear:
	/bin/rm -f decode_lz77 *.o
