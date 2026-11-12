# Compare GoSam with the reference values: make run
PROC ?= eepp
NLEG ?= 4
SRC ?= check_ref.f90
V=../virtual
FC=$(shell sh $(V)/config.sh -fortcom)
CFLAGS=$(shell sh $(V)/config.sh -cflags)
LDFLAGS=$(shell sh $(V)/config.sh -libs)

check_ref: $(SRC) momenta.inc
	sed -e 's/PROC_/$(PROC)_/g' -e 's/NLEG/$(NLEG)/' $(SRC) > check_ref_$(PROC).f90
	$(FC) $(CFLAGS) -cpp -o $@ check_ref_$(PROC).f90 $(LDFLAGS)

run: check_ref
	./check_ref

clean:
	rm -f check_ref check_ref_*.f90 *.o golem.log
