all: mgmres_prb

FC:=gfortran -std=f95

mgmres_prb: mgmres_prb.o mgmres.o
	$(FC) -o $@ $^

%.o: %.f90
	$(FC) -o $@ -c $^

.PHONY:: clean
clean:
	$(RM) mgmres_prb
	$(RM) mgmres_prb.o mgmres.o

.PHONY:: test
test: mgmres_prb mgmres_prb_output.txt
	./mgmres_prb | diff - mgmres_prb_output.txt
