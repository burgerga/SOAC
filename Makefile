#CC=pgf95
CC=gfortran

TRG1=smolarki
SRC1=smolarki.f95
OBJ1=$(SRC1:.f95=.out)

CLEANFLS=$(OBJ1) 

all: $(TRG1)

$(TRG1): $(OBJ1)

%.out: %.f95
	$(CC) $^ -o $@
	
clean:
	rm -f $(CLEANFLS) 

.PHONY: all clean
