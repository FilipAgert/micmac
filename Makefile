# Paths
DSRC = src
DOBJ = build
DEXE = app
DTEST = tests
DMOD = mod
DLIB = /lib/x86_64-linux-gnu
DSH = shell_correction
EXEN = main.exe
TEST_EXE = test.exe
FIT_EXEN = runfit.exe
FISS_EXEN = fiss_bar.exe
TAB_EXEN = write_tab.exe
HO_EXEN = ho.exe

# Flags
LIBS = -llapack -lblas -fopenmp
FLAGS = -Wall -O3 -I$(DOBJ) -ffree-line-length-none -fcheck=all -fbacktrace -g -fimplicit-none 
CC = gfortran $(FLAGS) -J$(DMOD) $(LIBS) -L$(DLIB) -c
CCL = gfortran -o

# Objects
OBJECTS = $(DOBJ)/constants.o $(DOBJ)/micmac.o $(DOBJ)/mass_table.o $(DOBJ)/fitting.o $(DOBJ)/optimise.o $(DOBJ)/def_ho.o $(DOBJ)/quadrule.o
TEST_OBJECTS = $(DOBJ)/test_micmac.o $(DOBJ)/test_utils.o $(DOBJ)/test_fitting.o $(DOBJ)/test_optimise.o $(DOBJ)/test_ho.o
MAIN_OBJ = $(DOBJ)/main.o
TEST_OBJ = $(DOBJ)/run_tests.o
FIT_OBJ = $(DOBJ)/runfit.o
FISS_OBJ = $(DOBJ)/run_fissbarr.o
HO_OBJ = $(DOBJ)/run_ho.o
TAB_OBJ = $(DOBJ)/table_writer.o
VPATH = $(DSRC):$(DTEST):$(DSRC)/$(DSH)

# Default target
all: main fit

$(DOBJ)/main.o: $(DSRC)/main.f90 $(DOBJ)/constants.o $(DOBJ)/micmac.o $(DOBJ)/mass_table.o #Which files does main depend on?
$(DOBJ)/fitting.o: $(DSRC)/fitting.f90 $(DOBJ)/constants.o $(DOBJ)/micmac.o $(DOBJ)/mass_table.o #Which files does fitting depend on?
$(DOBJ)/constants.o: $(DSRC)/constants.f90
$(DOBJ)/nucleus_module.o: $(DSRC)/nucleus_module.f90 $(DOBJ)/constants.o
$(DOBJ)/micmac.o: $(DSRC)/micmac.f90 $(DOBJ)/constants.o $(DOBJ)/optimise.o
$(DOBJ)/mass_table.o: $(DSRC)/mass_table.f90 $(DOBJ)/constants.o
$(DOBJ)/test_micmac.o: $(DOBJ)/micmac.o $(DOBJ)/constants.o $(DOBJ)/test_utils.o
$(DOBJ)/run_tests.o: $(DOBJ)/test_micmac.o $(DOBJ)/test_utils.o $(DOBJ)/test_fitting.o $(DOBJ)/test_optimise.o
$(DOBJ)/runfit.o: $(DOBJ)/constants.o $(DOBJ)/micmac.o $(DOBJ)/mass_table.o $(DOBJ)/fitting.o 
$(DOBJ)/test_fitting.o: $(DOBJ)/fitting.o
$(DOBJ)/run_fissbarr.o: $(DSRC)/run_fissbarr.f90 $(DOBJ)/constants.o $(DOBJ)/micmac.o $(DOBJ)/fitting.o
$(DOBJ)/table_writer.o: $(DSRC)/table_writer.f90 $(DOBJ)/constants.o $(DOBJ)/micmac.o
$(DOBJ)/test_ho.o: $(DTEST)/test_ho.f90 $(DOBJ)/def_ho.o $(DOBJ)/constants.o
$(DOBJ)/def_ho.o: $(DSRC)/$(DSH)/def_ho.f90 $(DOBJ)/constants.o $(DOBJ)/optimise.o $(DOBJ)/quadrule.o
$(DOBH)/run_ho.o: $(DSRC)/$(DSH)/run_ho.o $(DOBJ)/constants.o $(DOBJ)/def_ho.o
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

$(DEXE)/$(FISS_EXEN): $(FISS_OBJ) $(OBJECTS) | $(DEXE)
	$(CCL) $@ $(FISS_OBJ) $(OBJECTS) $(LIBS)

$(DEXE)/$(HO_EXEN): $(HO_OBJ) $(OBJECTS) | $(DEXE)
	$(CCL) $@ $(HO_OBJ) $(OBJECTS) $(LIBS)

$(DEXE)/$(TEST_EXE): $(TEST_OBJ) $(OBJECTS) $(TEST_OBJECTS) | $(DEXE)
	$(CCL) $@ $(TEST_OBJ) $(OBJECTS) $(TEST_OBJECTS) $(LIBS)


$(DEXE)/$(TAB_EXEN): $(TAB_OBJ) $(OBJECTS) | $(DEXE)
	$(CCL) $@ $(TAB_OBJ) $(OBJECTS) $(LIBS)

main: $(DEXE)/$(EXEN)

fit: $(DEXE)/$(FIT_EXEN)

fiss: $(DEXE)/$(FISS_EXEN)

buildtab: $(DEXE)/$(TAB_EXEN)

ho: $(DEXE)/$(HO_EXEN)

run: $(DEXE)/$(EXEN)
	$(DEXE)/$(EXEN)

test: $(DEXE)/$(TEST_EXE)
	$(DEXE)/$(TEST_EXE)

runfit: fit
	$(DEXE)/$(FIT_EXEN)

runfiss: fiss
	$(DEXE)/$(FISS_EXEN)

runtab: buildtab
	$(DEXE)/$(TAB_EXEN)

runho: ho
	$(DEXE)/$(HO_EXEN)

clean:
	rm -rf $(DOBJ)/*.o $(DEXE)/*.exe $(DMOD)/*.mod

.PHONY: clean run fit runfit main all
