all:		makeopac makeopac_distr makeopac_smoothed
makeopac:	make_scatmat.o bhmie.o Makefile
		gfortran make_scatmat.o bhmie.o -o makeopac
makeopac_distr:	make_scatmat_distr.o bhmie.o Makefile
		gfortran make_scatmat_distr.o bhmie.o -o makeopac_distr
makeopac_smoothed:	make_scatmat_smoothed.o bhmie.o Makefile
		gfortran make_scatmat_smoothed.o bhmie.o -o makeopac_smoothed
make_scatmat.o:	make_scatmat.f90 Makefile
		gfortran -c make_scatmat.f90
make_scatmat_distr.o:	make_scatmat_distr.f90 Makefile
		gfortran -c make_scatmat_distr.f90
make_scatmat_smoothed.o:	make_scatmat_smoothed.f90 Makefile
		gfortran -c make_scatmat_smoothed.f90
bhmie.o:	bhmie.f Makefile
		gfortran -c bhmie.f

clean:
	@rm -f	*.o *.mod *~
	@echo OBJECT and MODULE files removed.

cleanall:
	@rm -f	*.o *.mod *~ dustkapscatmat*.inp dustkappa*.inp makeopac makeopac_distr makeopac_smoothed *.out
