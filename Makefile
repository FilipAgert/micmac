# Paths
DSRC = src
DOBJ = build
DEXE = app
DTEST = tests
DMOD = mod
DLIB = /lib/x86_64-linux-gnu

EXEN = main.exe
TEST_EXE = test.exe
FIT_EXEN = runfit.exe

# Flags
LIBS = -llapack -lblas
FLAGS = -Wall -O3 -I$(DOBJ) -ffree-line-length-none -fcheck=all -fbacktrace -g -fimplicit-none
CC = gfortran $(FLAGS) -J$(DMOD) $(LIBS) -L$(DLIB) -c
CCL = gfortran -o

# Objects
OBJECTS = $(DOBJ)/constants.o $(DOBJ)/micmac.o $(DOBJ)/mass_table.o $(DOBJ)/table_writer.o $(DOBJ)/fitting.o
TEST_OBJECTS = $(DOBJ)/test_micmac.o $(DOBJ)/test_utils.o $(DOBJ)/test_fitting.o
MAIN_OBJ = $(DOBJ)/main.o
TEST_OBJ = $(DOBJ)/run_tests.o
FIT_OBJ = $(DOBJ)/runfit.o
VPATH = $(DSRC):$(DTEST)

# Default target
all: main fit

$(DOBJ)/main.o: $(DSRC)/main.f90 $(DOBJ)/constants.o $(DOBJ)/micmac.o $(DOBJ)/mass_table.o #Which files does main depend on?
$(DOBJ)/fitting.o: $(DSRC)/fitting.f90 $(DOBJ)/constants.o $(DOBJ)/micmac.o $(DOBJ)/mass_table.o $(DOBJ)/table_writer.o#Which files does fitting depend on?
$(DOBJ)/constants.o: $(DSRC)/constants.f90
$(DOBJ)/nucleus_module.o: $(DSRC)/nucleus_module.f90 $(DOBJ)/constants.o
$(DOBJ)/micmac.o: $(DSRC)/micmac.f90 $(DOBJ)/constants.o
$(DOBJ)/mass_table.o: $(DSRC)/mass_table.f90 $(DOBJ)/constants.o
$(DOBJ)/test_micmac.o: $(DOBJ)/micmac.o $(DOBJ)/constants.o $(DOBJ)/test_utils.o
$(DOBJ)/run_tests.o: $(DOBJ)/test_micmac.o $(DOBJ)/test_utils.o $(DOBJ)/test_fitting.o

# Ensure required directories exist
$(DOBJ) $(DEXE) $(DMOD) $(DTEST):
	mkdir -p $@

# Build rules
$(DOBJ)/%.o: %.f90 | $(DOBJ) $(DMOD)
	$(CC) $< -o $@


# Targets
$(DEXE)/$(EXEN): $(MAIN_OBJ) $(OBJECTS) | $(DEXE)
	$(CCL) $@ $(MAIN_OBJ) $(OBJECTS) $(LIBS)

$(DEXE)/$(FIT_EXEN): $(FIT_OBJ) $(OBJECTS) | $(DEXE)
	$(CCL) $@ $(FIT_OBJ) $(OBJECTS) $(LIBS)

$(DEXE)/$(TEST_EXE): $(TEST_OBJ) $(OBJECTS) $(TEST_OBJECTS) | $(DEXE)
	$(CCL) $@ $(TEST_OBJ) $(OBJECTS) $(TEST_OBJECTS) $(LIBS)

main: $(DEXE)/$(EXEN)

fit: $(DEXE)/$(FIT_EXEN)

run: $(DEXE)/$(EXEN)
	$(DEXE)/$(EXEN)

test: $(DEXE)/$(TEST_EXE)
	$(DEXE)/$(TEST_EXE)

runfit: fit
	$(DEXE)/$(FIT_EXEN)

clean:
	rm -rf $(DOBJ)/*.o $(DEXE)/*.exe $(DMOD)/*.mod

.PHONY: clean run fit runfit main all
