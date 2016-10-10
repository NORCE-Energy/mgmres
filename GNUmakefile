all: mgmres_prb

FC:=gfortran -std=f95

mgmres_prb: mgmres_prb.o mgmres.o
	$(FC) -o $@ $^

mgmres.o mgmres.mod: mgmres.f90
	$(FC) -o mgmres.o -c mgmres.f90

mgmres_prb.o: mgmres.mod mgmres_prb.f90
	$(FC) -o mgmres_prb.o -c mgmres_prb.f90

.PHONY:: clean
clean:
	$(RM) mgmres_prb
	$(RM) mgmres_prb.o mgmres.o mgmres.mod

.PHONY:: test
test: mgmres_prb mgmres_prb_output.txt
	./mgmres_prb | diff - mgmres_prb_output.txt
