# Compare GoSam with the reference values: make run
# without vacuum polarisation (process in ../virtual_novp): make run NOVP=1
# (python3 check_vp.py runs this and adds the vacuum-polarisation subtraction)
PROC ?= eeppA
NLEG ?= 5
SRC ?= ../../eepp_sqed/check_ref/check_ref.f90
ifdef NOVP
V=../virtual_novp
EXE=check_ref_novp
DEFS=-DNOVP
else
V=../virtual
EXE=check_ref
DEFS=
endif
FC=$(shell sh $(V)/config.sh -fortcom)
CFLAGS=$(shell sh $(V)/config.sh -cflags)
LDFLAGS=$(shell sh $(V)/config.sh -libs)

$(EXE): $(SRC) momenta.inc FORCE
	sed -e 's/PROC_/$(PROC)_/g' -e 's/NLEG/$(NLEG)/' $(SRC) > $(EXE)_$(PROC).f90
	$(FC) $(CFLAGS) -cpp $(DEFS) -o $@ $(EXE)_$(PROC).f90 $(LDFLAGS)

run: $(EXE)
	./$(EXE)

.PHONY: run clean FORCE
FORCE:

clean:
	rm -rf check_ref check_ref_novp check_ref*_*.f90 *.dSYM *.o golem.log
