#CC=pgf95
CC=gfortran
CFLAGS=-g

TRG1=smolarki
SRC1=smolarki.f95
OBJ1=$(SRC1:.f95=.out)

TRG2=smolarki2d
SRC2=smolarki2d.f95
OBJ2=$(SRC2:.f95=.out)

TRG3=smolarki3d
SRC3=smolarki3d.f95
OBJ3=$(SRC3:.f95=.out)

CLEANFLS=$(OBJ1) $(OBJ2) $(OBJ3) $(TRG1) $(TRG2) $(TRG3)
#CLEANFLS=$(OBJ1) $(TGR1)

all: $(TRG1) $(TRG2) $(TRG3)
#all: $(TRG1)

$(TRG1): $(OBJ1)
$(TRG2): $(OBJ2)
$(TRG3): $(OBJ3)

%.out: %.f95
	$(CC) $(CFLAGS) $^ -o $@

.PHONY: all clean
