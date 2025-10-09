# --- Makefile simple pour gfortran ---
FC      = gfortran
FFLAGS  = -O2 -Wall -std=f2008
TARGET  = run
OBJS    = precision.o donnees.o schema_rusanov.o main.o

all: $(TARGET)

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS)

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: all run clean
