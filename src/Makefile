# Makefile for Radex
# Uncomment the FC/FFLAGS combination of your choice
# Option 1:
#FC      = g77
#FFLAGS += -fno-automatic -Wall -O2
# Option 2:
#FC     = ifort
#f2pyFC = intelem #this must match FC (ifort=intelem, gfortran=gnu95)
#FFLAGS += -g -fPIC -fmax-errors=3 -fcheck=all -fbacktrace
#FFLAGS += -g -fPIC

# Option 3:
FC      = gfortran
#FFLAGS += -O2
FFLAGS += -g -fPIC
# Option 4:
#FC = g95
#FFLAGS += -Wall -O2
# there may be other options of course ...

SRCS = types slatec commondata solver io background main
OBJS = $(SRCS:=.o)
EXEC    = nradex
BINDIR  = ../bin/

all: $(EXEC) install clean

$(EXEC): $(OBJS)
	$(FC) $(FFLAGS) -o $(EXEC) $(OBJS)

.SUFFIXES: .o .f90

.f90.o:
	$(FC) -c $(FCFLAGS) $<

install: $(EXEC)
	install -m 755 -p -s $(EXEC) $(BINDIR)

clean:
	- rm *.o *.mod
	- rm $(EXEC)
