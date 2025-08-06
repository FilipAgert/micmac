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
PAIR_EXEN = run_pair.exe
FISS_EXEN = fiss_bar.exe
TAB_EXEN = write_tab.exe
SH_EXEN = shell
ST_EXEN = strut

# Flags
LIBS = -llapack -lblas -fopenmp
FLAGS = -O2 -I$(DOBJ) -ffree-line-length-none -fcheck=all -fbacktrace -g -fimplicit-none -fno-omit-frame-pointer
CC = gfortran $(FLAGS) -J$(DMOD) $(LIBS) -L$(DLIB) -c
CCL = gfortran -o

# Objects
OBJECTS = $(DOBJ)/constants.o $(DOBJ)/micmac.o $(DOBJ)/mass_table.o $(DOBJ)/fitting.o $(DOBJ)/optimise.o $(DOBJ)/def_ho.o $(DOBJ)/quad.o $(DOBJ)/hamiltonian.o $(DOBJ)/strutinsky.o $(DOBJ)/brent.o $(DOBJ)/pairing.o
TEST_OBJECTS = $(DOBJ)/test_micmac.o $(DOBJ)/test_utils.o $(DOBJ)/test_fitting.o $(DOBJ)/test_optimise.o $(DOBJ)/test_ho.o
MAIN_OBJ = $(DOBJ)/main.o
TEST_OBJ = $(DOBJ)/run_tests.o
FIT_OBJ = $(DOBJ)/runfit.o
FISS_OBJ = $(DOBJ)/run_fissbarr.o
SH_OBJ = $(DOBJ)/run_shell.o
TAB_OBJ = $(DOBJ)/table_writer.o
PAIR_OBJ = $(DOBJ)/run_pair.o
ST_OBJ =$(DOBJ)/run_strut.o
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
$(DOBJ)/test_fitting.o: $(DOBJ)/fitting.o $(DOBJ)/constants.o
$(DOBJ)/run_fissbarr.o: $(DSRC)/run_fissbarr.f90 $(DOBJ)/constants.o $(DOBJ)/micmac.o $(DOBJ)/fitting.o $(DOBJ)/constants.o
$(DOBJ)/table_writer.o: $(DSRC)/table_writer.f90 $(DOBJ)/constants.o $(DOBJ)/micmac.o 
$(DOBJ)/test_ho.o: $(DTEST)/test_ho.f90 $(DOBJ)/def_ho.o $(DOBJ)/constants.o
$(DOBJ)/def_ho.o: $(DSRC)/$(DSH)/def_ho.f90 $(DOBJ)/constants.o $(DOBJ)/optimise.o $(DOBJ)/quad.o $(DOBJ)/micmac.o
$(DOBJ)/run_ho.o: $(DSRC)/$(DSH)/run_ho.f90 $(DOBJ)/constants.o $(DOBJ)/def_ho.o $(DOBJ)/hamiltonian.o
$(DOBJ)/run_strut.o: $(DSRC)/$(DSH)/run_strut.f90 $(DOBJ)/constants.o $(DOBJ)/def_ho.o $(DOBJ)/hamiltonian.o $(DOBJ)/strutinsky.o $(DOBJ)/brent.o
$(DOBJ)/strutinsky.o: $(DSRC)/$(DSH)/strutinsky.f90 $(DOBJ)/constants.o $(DOBJ)/def_ho.o $(DOBJ)/hamiltonian.o $(DOBJ)/brent.o $(DOBJ)/pairing.o
$(DOBJ)/hamiltonian.o: $(DSRC)/$(DSH)/hamiltonian.f90 $(DOBJ)/constants.o $(DOBJ)/def_ho.o $(DOBJ)/optimise.o $(DOBJ)/quad.o
$(DOBJ)/pairing.o: $(DSRC)/$(DSH)/pairing.f90 $(DOBJ)/constants.o
$(DOBJ)/run_pair.o: $(DSRC)/$(DSH)/run_pair.f90 $(DOBJ)/pairing.o $(DOBJ)/constants.o
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

$(DEXE)/$(SH_EXEN): $(SH_OBJ) $(OBJECTS) | $(DEXE)
	$(CCL) $@ $(SH_OBJ) $(OBJECTS) $(LIBS)

$(DEXE)/$(ST_EXEN): $(ST_OBJ) $(OBJECTS) | $(DEXE)
	$(CCL) $@ $(ST_OBJ) $(OBJECTS) $(LIBS)

$(DEXE)/$(TEST_EXE): $(TEST_OBJ) $(OBJECTS) $(TEST_OBJECTS) | $(DEXE)
	$(CCL) $@ $(TEST_OBJ) $(OBJECTS) $(TEST_OBJECTS) $(LIBS)


$(DEXE)/$(TAB_EXEN): $(TAB_OBJ) $(OBJECTS) | $(DEXE)
	$(CCL) $@ $(TAB_OBJ) $(OBJECTS) $(LIBS)

$(DEXE)/$(PAIR_EXEN): $(PAIR_OBJ) $(OBJECTS) | $(DEXE)
	$(CCL) $@ $(PAIR_OBJ) $(OBJECTS) $(LIBS)

main: $(DEXE)/$(EXEN)

fit: $(DEXE)/$(FIT_EXEN)

fiss: $(DEXE)/$(FISS_EXEN)

buildtab: $(DEXE)/$(TAB_EXEN)

sh: $(DEXE)/$(SH_EXEN)

st: $(DEXE)/$(ST_EXEN)

pair: $(DEXE)/$(PAIR_EXEN)


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

runsh: sh
	$(DEXE)/$(SH_EXEN)

runst: st
	$(DEXE)/$(ST_EXEN)

runpair: pair
	$(DEXE)/$(PAIR_EXEN)

runcallgrind:
	@name="callgrind.out.$(SH_EXEN)"; \
	n=0; \
	filename="$$name"; \
	while [ -e "$$filename" ]; do \
		n=$$((n+1)); \
		filename="$$name\_$$n"; \
	done; \
	echo "Saving callgrind output to $$filename"; \
	valgrind --tool=callgrind --callgrind-out-file=$$filename $(DEXE)/$(SH_EXEN); \
	kcachegrind $$filename

clean:
	rm -rf $(DOBJ)/*.o $(DEXE)/*.exe $(DMOD)/*.mod

.PHONY: clean run fit runfit main all
