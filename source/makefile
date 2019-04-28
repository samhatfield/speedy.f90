# Fortran compiler
FC=gfortran

# Set variables depending on target
default : COMPOPTS = $(OPT) $(BASE)
profile : COMPOPTS = -pg $(BASE)

# Default
default : base_target

# For profiling
profile : base_target

# Base compiler options (always used)
BASE=-fdefault-real-8 -fconvert=swap -Wall

# Optimisation flags (disabled for debugging, profiling etc.)
OPT=-Ofast

# # Location of NetCDF module (netcdf.mod)
INC=-I$(NETCDF)

# Library flags
LIB=-lnetcdff -lnetcdf

FILES= \
	   boundaries.o \
       coupler.o \
	   date.o \
	   diagnostics.o \
	   dynamical_constants.o \
	   forcing.o \
	   fourier.o \
	   geometry.o \
 	   geopotential.o \
	   horizontal_diffusion.o \
	   humidity.o \
	   implicit.o \
       initialization.o \
	   land_model.o \
	   legendre.o \
	   longwave_radiation.o \
	   matrix_inversion.o \
	   mod_cpl_bcinterp.o \
	   mod_dyncon2.o \
	   input_output.o \
	   mod_radcon.o \
	   params.o \
	   physics.o \
	   physical_constants.o \
       ppo_dmflux.o \
	   precipitation.o \
	   prognostics.o \
	   sea_model.o \
	   shortwave_radiation.o \
	   spectral.o \
       fftpack.o \
	   sppt.o \
	   surface_fluxes.o \
	   tendencies.o \
	   time_stepping.o \
	   vertical_diffusion.o

%.o: %.f90
	$(FC) $(COMPOPTS) -c $< $(INC)

base_target: $(FILES) speedy.o
	$(FC) $(COMPOPTS) $(FILES) speedy.o -o speedy $(LIB)

.PHONY: clean
clean:
	rm -f *.o *.mod

speedy.o               : params.o date.o input_output.o shortwave_radiation.o time_stepping.o\
                         diagnostics.o
boundaries.o           : physical_constants.o params.o input_output.o spectral.o
coupler.o              : land_model.o sea_model.o
fourier.o              : params.o geometry.o fftpack.o
geometry.o             : params.o physical_constants.o
geopotential.o         : params.o physical_constants.o geometry.o
horizontal_diffusion.o : params.o dynamical_constants.o
humidity.o             : params.o
implicit.o             : params.o dynamical_constants.o physical_constants.o mod_dyncon2.o\
                         geometry.o horizontal_diffusion.o matrix_inversion.o
initialization.o       : coupler.o params.o date.o input_output.o time_stepping.o boundaries.o\
                         spectral.o sea_model.o physics.o geopotential.o prognostics.o forcing.o
forcing.o              : dynamical_constants.o shortwave_radiation.o params.o \
                         physical_constants.o boundaries.o date.o land_model.o mod_radcon.o\
						 surface_fluxes.o date.o sea_model.o longwave_radiation.o humidity.o\
						 horizontal_diffusion.o
land_model.o           : params.o date.o mod_cpl_bcinterp.o input_output.o boundaries.o
legendre.o             : params.o physical_constants.o geometry.o
diagnostics.o          : params.o spectral.o
mod_dyncon2.o          : params.o
prognostics.o          : params.o dynamical_constants.o physical_constants.o geometry.o\
                         boundaries.o diagnostics.o spectral.o input_output.o
input_output.o         : params.o  physical_constants.o date.o spectral.o geometry.o
physical_constants.o   : params.o
mod_radcon.o           : params.o
physics.o              : params.o coupler.o physical_constants.o boundaries.o land_model.o\
                         sea_model.o sppt.o precipitation.o surface_fluxes.o vertical_diffusion.o\
						 shortwave_radiation.o longwave_radiation.o humidity.o geometry.o
longwave_radiation.o   : params.o physical_constants.o mod_radcon.o geometry.o
sea_model.o            : params.o input_output.o boundaries.o geometry.o mod_cpl_bcinterp.o\
 						 date.o
shortwave_radiation.o  : params.o mod_radcon.o geometry.o
surface_fluxes.o       : params.o physical_constants.o mod_radcon.o land_model.o humidity.o
ppo_dmflux.o           : params.o physical_constants.o land_model.o sea_model.o physics.o\
                         mod_radcon.o
precipitation.o        : params.o physical_constants.o
spectral.o             : params.o physical_constants.o legendre.o fourier.o
sppt.o                 : params.o physical_constants.o spectral.o
tendencies.o           : params.o implicit.o prognostics.o physical_constants.o geometry.o\
                         mod_dyncon2.o physics.o spectral.o geopotential.o
time_stepping.o        : dynamical_constants.o params.o prognostics.o tendencies.o\
                         horizontal_diffusion.o
vertical_diffusion     : params.o physical_constants.o geometry.o
