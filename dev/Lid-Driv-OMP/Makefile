
# Inspired on $UWHPSC/codes/fortran/multifile1/Makefile{5,7}

FC = gfortran    
FFLAGS = -O3
LFLAGS =
SOURCES = $(wildcard *.f90)
SOURCES = lidDrivOMP.f90
OBJECTS = $(subst .f90,.o,$(SOURCES))

.PHONY: test clean help

lidDrivOMP.out: lidDrivOMP.exe
	./lidDrivOMP.exe

lidDrivOMP.exe: $(OBJECTS)
	$(FC) $(LFLAGS) $(OBJECTS) -o lidDrivOMP.exe

%.o : %.f90
	$(FC) $(FFLAGS) -c $< 

clean:
	rm -f $(OBJECTS) lidDrivOMP.exe

test:
	@echo "Sources are: " $(SOURCES)
	@echo "Objects are: " $(OBJECTS)

help:
	@echo "Valid targets:      "
	@echo "  lidDrivOMP.out    "
	@echo "  lidDrivOMP.exe    "
	@echo "  lidDrivOMP.o      "
	@echo "  clean, test, help "

