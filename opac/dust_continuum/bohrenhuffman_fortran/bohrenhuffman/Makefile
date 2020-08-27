makeopac:	make_ca_cs_g.o bhmie.o Makefile
		gfortran make_ca_cs_g.o bhmie.o -o makeopac
make_ca_cs_g.o:	make_ca_cs_g.f90 Makefile
		gfortran -c make_ca_cs_g.f90
bhmie.o:	bhmie.f Makefile
		gfortran -c bhmie.f

clean:
	@rm -f	*.o *.mod *~
	@echo OBJECT and MODULE files removed.

cleanall:
	@rm -f	*.o *.mod *~ dustkappa*.inp makeopac
