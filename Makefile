CC=mpicc
CFLAGS= 
CFLAGS= -I$(HOME)/local/include 
#LDFLAGS=-lm
LDFLAGS=-L$(HOME)/local/lib -lm -lmpi -lgsl -lgslcblas #-lfftw3 
EXEC=swot
SRC=  main.c
OBJ= $(SRC:.c=.o)

all: $(EXEC)

$(EXEC) : $(OBJ)
	$(CC) -o $@ $^ $(LDFLAGS)

main.o: main.h

%.o: %.c
	$(CC) -o $@ -c $< $(CFLAGS)

.PHONY: clean mrproper

clean:
	rm -rf *.o $(EXEC).tar

mrproper: clean
	rm -rf $(EXEC)

tar:
	tar cvf $(EXEC).tar Makefile main.c main.h README $(EXEC)

gzip:
	gzip $(EXEC).tar