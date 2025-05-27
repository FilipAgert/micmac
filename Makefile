# Paths
DSRC = src
DOBJ = build
DEXE = app
DTEST = tests
DMOD = mod

EXEN = main.exe
FIT_EXEN = runfit.exe

# Flags
FLAGS = -Wall -O3 -I$(DMOD) -ffree-line-length-none -fcheck=all -fbacktrace -g
CC = gfortran $(FLAGS) -J$(DMOD) -c
CCL = gfortran -o

# Objects
OBJECTS = $(DOBJ)/constants.o $(DOBJ)/nucleus_module.o $(DOBJ)/micmac.o $(DOBJ)/mass_table.o 
MAIN_OBJ = $(DOBJ)/main.o
FIT_OBJ = $(DOBJ)/fitting.o

# Default target
all: main fit

$(DOBJ)/main.o: $(DSRC)/main.f90 $(DOBJ)/constants.o $(DOBJ)/nucleus_module.o $(DOBJ)/micmac.o $(DOBJ)/mass_table.o #Which files does main depend on?
$(DOBJ)/fitting.o: $(DSRC)/fitting.f90 $(DOBJ)/constants.o $(DOBJ)/micmac.o $(DOBJ)/mass_table.o#Which files does fitting depend on?
$(DOBJ)/constants.o: $(DSRC)/constants.f90
$(DOBJ)/nucleus_module.o: $(DSRC)/nucleus_module.f90 $(DOBJ)/constants.o
$(DOBJ)/micmac.o: $(DSRC)/micmac.f90 $(DOBJ)/constants.o
$(DOBJ)/mass_table.o: $(DSRC)/mass_table.f90 $(DOBJ)/constants.o

# Ensure required directories exist
$(DOBJ) $(DEXE) $(DMOD):
	mkdir -p $@

# Build rules
$(DOBJ)/%.o: $(DSRC)/%.f90 | $(DOBJ) $(DMOD)
	$(CC) $< -o $@

# Targets
$(DEXE)/$(EXEN): $(MAIN_OBJ) $(OBJECTS) | $(DEXE)
	$(CCL) $@ $(MAIN_OBJ) $(OBJECTS)

$(DEXE)/$(FIT_EXEN): $(FIT_OBJ) $(OBJECTS) | $(DEXE)
	$(CCL) $@ $(FIT_OBJ) $(OBJECTS)

main: $(DEXE)/$(EXEN)

fit: $(DEXE)/$(FIT_EXEN)

run:
	$(DEXE)/$(EXEN)

runfit: fit
	$(DEXE)/$(FIT_EXEN)

clean:
	rm -rf $(DOBJ)/*.o $(DEXE)/*.exe $(DMOD)/*.mod

.PHONY: clean run fit runfit main all
