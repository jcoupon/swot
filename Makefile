# Makefile for swot

# Where cfitsio library is installed
CFITSIO =  #/usr/local


# Where GSL library is installed
GSL =  #/usr/local

# Where MPI is installed
MPI = #/opt/openmpi-1.8.6_clang

# compiler options
ifneq ($(MPI), )
	MPICC       = $(MPI)/bin/mpicc
else
	MPICC       = mpicc
endif

CFLAGS      = -Iinclude # -use-asm # for old icc versions
LDFLAGS     =  -lgsl -lgslcblas -lm -lcfitsio
MPI_CFLAGS  =
MPI_LFLAGS  =
RM          = rm -f
EXEC        = bin/swot
SRC         = main.c
OBJ         = $(SRC:.c=.o)

# source files
SRCS    = utils.c tree.c init.c  correlators.c main.c
OBJS    = $(SRCS:.c=.o)

# Headers for libraries
ifneq ($(CFITSIO), )
	CFLAGS     +=  -I$(CFITSIO)/include
	LDFLAGS    +=  -L$(CFITSIO)/lib
endif

ifneq ($(GSL), )
	CFLAGS     +=  -I$(GSL)/include
	LDFLAGS    +=  -L$(GSL)/lib
endif

ifneq ($(MPI), )
	MPI_CFLAGS +=  -I$(MPI)/include
	MPI_LFLAGS +=  -L$(MPI)/lib
endif

.PHONY: all
all: $(EXEC)

vpath %.h include
vpath %.c src

$(EXEC):  $(OBJS)
	$(MPICC)  -o $@ $^ $(CFLAGS) $(LDFLAGS) $(MPI_CFLAGS) $(MPI_LFLAGS)

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
