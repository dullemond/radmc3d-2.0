all:		makeopac_smoothed
makeopac_smoothed:	make_scatmat_smoothed.o bhmie.o Makefile
		gfortran make_scatmat_smoothed.o bhmie.o -o makeopac_smoothed
make_scatmat_smoothed.o:	make_scatmat_smoothed.f90 Makefile
		gfortran -c make_scatmat_smoothed.f90
bhmie.o:	bhmie.f Makefile
		gfortran -c bhmie.f

clean:
	@rm -f	*.o *.mod *~
	@echo OBJECT and MODULE files removed.

cleanall:
	@rm -f	*.o *.mod *~ *.inp makeopac_smoothed 

cleanopac:
	@rm -f	dustkapscatmat*.inp dustkappa*.inp
