# Makefile for swot

# compiler options
MPICC       = mpicc
CFLAGS      = # -use-asm # for old icc versions
LDFLAGS      = 
MPI_CFLAGS  = 
MPI_LFLAGS  = 
RM          = rm -f
EXEC        = bin/swot
SRC         = main.c
OBJ         = $(SRC:.c=.o)

# Where GSL library is installed
GSL = /usr/local

# Where MPI is installed
MPI = # /opt/openmpi-1.8.5/

# source files
SRCS    = main.c
OBJS    = $(SRCS:.c=.o)

# Headers for libraries
CFLAGS     +=  -Iinclude  -I$(GSL)/include
LDFLAGS     += -lm -lgsl -lgslcblas  -L$(GSL)/lib 
MPI_CFLAGS +=  -I$(MPI)/include
MPI_LFLAGS +=  -L$(MPI)/lib 

.PHONY: all
all: $(EXEC)

vpath %.h include
vpath %.c src

$(EXEC):  $(OBJS)
	$(MPICC) $(CFLAGS) $(LDFLAGS) $(MPI_CFLAGS) $(MPI_LFLAGS) -o $@ $^

%.o:  %.c 
	$(MPICC) -c -o $@ $< $(CFLAGS) $(MPI_CFLAGS)

%.h:

.PHONY: clean
clean:
	-${RM} ${OBJS}

#tar:
#	tar cvf $(EXEC).tar Makefile main.c main.h README $(EXEC)

#gzip:
#	gzip $(EXEC).tar
