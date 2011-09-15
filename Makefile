#CC=pgf95
CC=gfortran

TRG1=exc2_4
SRC1=exc2_4.f90
OBJ1=$(SRC1:.f90=.out)

TRG2=exc3_2
SRC2=exc3_2.f90
OBJ2=$(SRC2:.f90=.out)

TRG3=exc3_4
SRC3=exc3_4.f90
OBJ3=$(SRC3:.f90=.out)

CLEANFLS=$(OBJ1) $(OBJ2) $(OBJ3)

all: $(TRG1) $(TRG2) $(TRG3)

$(TRG1): $(OBJ1)
$(TRG2): $(OBJ2)
$(TRG3): $(OBJ3)

%.out: %.f90
	$(CC) $^ -o $@
	
clean:
	rm -f $(CLEANFLS) 

.PHONY: all clean
