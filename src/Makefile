SRC = .
#SRC = /Users/cornelisdullemond/science/prog/radmc-3d/version_0.41/src

# Optimize
OPTIM = -O2 -fopenmp
# OPTIM = -O2 

# The system defines used in executing make
# FF = /Users/dullemon/libraries/g95-install/bin/powerpc-apple-darwin6.8-g95

# For using the OpenMP parallel version of RADMC-3D one has to add '-fopenmp' in
# the next command line (see also the RADMC-3D manual).
FF = gfortran
#FF = g95

# Boundary check
# BCHECK   = -fbounds-check

# Glueing all the things together
FFLAGS   = $(OPTIM) $(BCHECK)
#FFLAGS   = $(OPTIM)
#-Wall
LFLAGS   = $(OPTIM) 
#-Wall

#
# Names
#
SRCMAIN        = $(SRC)/main.f90
SRCRTGLOBAL    = $(SRC)/rtglobal_module.f90
SRCMONTE       = $(SRC)/montecarlo_module.f90
SRCDUST        = $(SRC)/dust_module.f90
SRCPOL         = $(SRC)/polarization_module.f90
SRCQUANTUM     = $(SRC)/quantum_module.f90
SRCMATH        = $(SRC)/mathroutines_module.f90
SRCIOPUT       = $(SRC)/ioput_module.f90
SRCSTARS       = $(SRC)/stars_module.f90
SRCAMR         = $(SRC)/amr_module.f90
SRCAMRRAY      = $(SRC)/amrray_module.f90
SRCCONSTANTS   = $(SRC)/constants_module.f90
SRCCAMERA      = $(SRC)/camera_module.f90
SRCLINES       = $(SRC)/lines_module.f90
SRCGASCONT     = $(SRC)/gascontinuum_module.f90
SRCNAMELIST    = $(SRC)/namelist_module.f90
SRCSOURCES     = $(SRC)/sources_module.f90
SRCUSERDEF     = ./userdef_module.f90


OBJ            = main.o rtglobal_module.o montecarlo_module.o dust_module.o \
	         quantum_module.o mathroutines_module.o ioput_module.o      \
                 stars_module.o amr_module.o amrray_module.o                \
	         constants_module.o camera_module.o lines_module.o          \
	         namelist_module.o userdef_module.o gascontinuum_module.o   \
	         polarization_module.o sources_module.o
GMODULES       = rtglobal_module.o mathroutines_module.o amr_module.o       \
	         namelist_module.o constants_module.o
HMODULES       = $(GMODULES) polarization_module.o
HIMODULES      = $(HMODULES) amrray_module.o dust_module.o 
IMODULES       = $(HIMODULES) stars_module.o
DMODULES       = $(IMODULES) quantum_module.o lines_module.o ioput_module.o \
	                     gascontinuum_module.o
UMODULES       = $(DMODULES) montecarlo_module.o 
SMODULES       = $(UMODULES) userdef_module.o
TMODULES       = $(SMODULES) sources_module.o
MODULES        = $(TMODULES) camera_module.o

#################################################
#                   RULES                       #
#################################################


all:	            radmc3d

radmc3d:            $(OBJ) Makefile
	            $(FF) $(LFLAGS) $(OBJ) $(LIBS) -o $@ 

main.o:             $(SRCMAIN) $(MODULES) Makefile
	            $(FF) -c $(FFLAGS)  $(SRCMAIN) -o $@

rtglobal_module.o:  $(SRCRTGLOBAL) amr_module.o Makefile
	            $(FF) -c $(FFLAGS)  $(SRCRTGLOBAL) -o $@

montecarlo_module.o: $(SRCMONTE) $(DMODULES) Makefile
	            $(FF) -c $(FFLAGS)  $(SRCMONTE) -o $@

dust_module.o:      $(SRCDUST) $(HMODULES) Makefile
	            $(FF) -c $(FFLAGS)  $(SRCDUST) -o $@

polarization_module.o: $(SRCPOL) $(GMODULES) Makefile
	            $(FF) -c $(FFLAGS)  $(SRCPOL) -o $@

quantum_module.o:   $(SRCQUANTUM) $(GMODULES) Makefile
	            $(FF) -c $(FFLAGS)  $(SRCQUANTUM) -o $@

mathroutines_module.o: $(SRCMATH) Makefile
	            $(FF) -c $(FFLAGS)  $(SRCMATH) -o $@

constants_module.o: $(SRCCONSTANTS) Makefile
	            $(FF) -c $(FFLAGS)  $(SRCCONSTANTS) -o $@

stars_module.o:     $(SRCSTARS) $(HIMODULES) Makefile
	            $(FF) -c $(FFLAGS)  $(SRCSTARS) -o $@

ioput_module.o:     $(SRCIOPUT) $(IMODULES) Makefile
	            $(FF) -c $(FFLAGS)  $(SRCIOPUT) -o $@

amr_module.o:       $(SRCAMR) constants_module.o Makefile
	            $(FF) -c $(FFLAGS)  $(SRCAMR) -o $@

amrray_module.o:    $(SRCAMRRAY) constants_module.o Makefile
	            $(FF) -c $(FFLAGS)  $(SRCAMRRAY) -o $@

camera_module.o:    $(SRCCAMERA) $(SMODULES) Makefile
	            $(FF) -c $(FFLAGS)  $(SRCCAMERA) -o $@

lines_module.o:     $(SRCLINES) $(IMODULES) ioput_module.o dust_module.o quantum_module.o  Makefile
	            $(FF) -c $(FFLAGS)  $(SRCLINES) -o $@

gascontinuum_module.o:  $(SRCGASCONT) $(IMODULES) ioput_module.o Makefile
	            $(FF) -c $(FFLAGS)  $(SRCGASCONT) -o $@

sources_module.o:   $(SRCSOURCES) $(DMODULES) Makefile
	            $(FF) -c $(FFLAGS)  $(SRCSOURCES) -o $@

namelist_module.o:  $(SRCNAMELIST) Makefile
	            $(FF) -c $(FFLAGS)  $(SRCNAMELIST) -o $@

userdef_module.o:   $(SRCUSERDEF) $(UMODULES) Makefile
	            $(FF) -c $(FFLAGS)  $(SRCUSERDEF) -o $@


install:	    radmc3d
	            ./install.perl

clean:
	@rm -f	*.o *.mod
	@echo OBJECT and MODULE files removed.

cleanmodel:
	@rm -f	*.o *.mod *.dat *.out *.info *.used                 \
	   wavelength_micron.inp dustopac.inp dust_density.inp      \
	   dust_density.uinp *.uout *.udat numberdens_*inp          \
	   amr_grid.inp stars.inp radmc3d.inp gas_temperature.*inp  \
	   gas_velocity.inp lines.inp microturbulence.*             \
	   *.binp *.bout *.bdat *.vtk                               \
	   electrons_numdens.*inp ions_numdens.*inp
	@echo All .o, .mod, .out files and most .inp files removed.

cleanall:
	@rm -f	*.o *.mod *.dat *.out *.info *.used                 \
	   wavelength_micron.inp dustopac.inp dust_density.inp      \
	   dust_density.uinp *.uout *.udat numberdens_*inp          \
	   amr_grid.inp stars.inp radmc3d.inp image_script.pro      \
	   radmc3d Makefile~ *.pro~ README*~ *.f90~ *.inp~          \
	   gas_velocity.inp lines.inp  microturbulence.*            \
	   gas_temperature.*inp electrons_numdens.*inp              \
	   *.binp *.bout *.bdat *.vtk                               \
	   ions_numdens.*inp
	@echo Directory cleaned to basic
