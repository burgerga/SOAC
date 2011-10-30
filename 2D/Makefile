#CC=pgf95
CC=gfortran-4.5
CFLAGS=-g

TRG2=smolarki2d
SRC2=smolarki2d.f95
OBJ2=$(SRC2:.f95=.out)

CLEANFLS= $(OBJ2) 

all: $(TRG2)

$(TRG2): $(OBJ2)

%.out: %.f95
	$(CC) $(CFLAGS) $^ -o $@

.PHONY: all clean
