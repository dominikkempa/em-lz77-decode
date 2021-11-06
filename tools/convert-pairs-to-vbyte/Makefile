SHELL = /bin/sh
CC = g++
CFLAGS = -Wall -Wextra -pedantic -Wshadow -funroll-loops -std=c++0x -O3 -DNDEBUG -march=native
#CFLAGS = -Wall -Wextra -pedantic -Wshadow -pthread -std=c++0x -g2

all: convert
convert:
	$(CC) $(CFLAGS) -o convert utils.cpp main.cpp
clean:
	/bin/rm -f *.o
nuclear:
	/bin/rm -f convert *.o
