#CC=pgf95
CC=gfortran
CFLAGS=-g

TRG1=smolarki
SRC1=smolarki.f95
OBJ1=$(SRC1:.f95=.out)

TRG2=smolarki2d
SRC2=smolarki2d.f95
OBJ2=$(SRC2:.f95=.out)

CLEANFLS=$(OBJ1) $(OBJ2) $(TRG1) $(TRG2)

all: $(TRG1) $(TRG2)

$(TRG1): $(OBJ1)
$(TRG2): $(OBJ2)

%.out: %.f95
	$(CC) $(CFLAGS) $^ -o $@

<<<<<<< HEAD
%.out: %.f95
	$(CC) $^ -o $@
	
clean:
	rm -f $(CLEANFLS) 
=======
clean: rm -f $(CLEANFLS)
>>>>>>> 7740b6860ffee05f5f5494f5cf1ddb06b2e4b67e

.PHONY: all clean
