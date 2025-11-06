# --- Makefile simple pour gfortran ---
FC      = gfortran
FFLAGS  = -O2 -Wall -std=f2008
TARGET  = run
OBJS    = precision.o fonctions.o donnees.o schema_rusanov.o main.o

all: $(TARGET)

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS)



clean:
	rm -f $(OBJS) $(TARGET) *.dat *.mod *.o

.PHONY: all run clean
