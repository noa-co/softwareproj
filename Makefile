CC=gcc
CFLAGS= -ansi -Wall -Wextra -Werror -pedantic-errors -lm
OBJ = spkmeans

spkmeans: spkmeans.o sp_io_utils.o sp_kmeansalgo.o
	$(CC) -o spkmeans spkmeans.o sp_io_utils.o sp_kmeansalgo.o $(CFLAGS)

spkmeans.o: spkmeans.c
	$(CC) -c spkmeans.c spkmeans.h $(CFLAGS)

sp_io_utils.o: sp_io_utils.c
	$(CC) -c sp_io_utils.c spkmeans.h $(CFLAGS)

sp_kmeansalgo.o: sp_kmeansalgo.c
	$(CC) -c sp_kmeansalgo.c spkmeans.h $(CFLAGS)