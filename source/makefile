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
BASE=-fdefault-real-8 -fconvert=swap

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
       dyn_geop.o \
       dyn_grtend.o \
       dyn_sptend.o \
	   humidity.o \
	   implicit.o \
       initialization.o \
       ini_indyns.o \
       ini_invars.o \
       ini_fordate.o \
	   land_model.o \
	   longwave_radiation.o \
	   mod_atparam.o \
	   mod_cpl_bcinterp.o \
	   mod_dyncon0.o \
	   mod_dyncon1.o \
	   mod_dyncon2.o \
	   mod_hdifcon.o \
	   input_output.o \
	   mod_radcon.o \
	   mod_tsteps.o \
	   physics.o \
	   physical_constants.o \
       ppo_dmflux.o \
	   precipitation.o \
	   prognostics.o \
	   sea_model.o \
	   shortwave_radiation.o \
	   spectral.o \
       spe_matinv.o \
       spe_spectral.o \
       spe_subfft_fftpack.o \
	   sppt.o \
	   surface_fluxes.o \
	   time_stepping.o \
	   vertical_diffusion.o

%.o: %.f90
	$(FC) $(COMPOPTS) -c $< $(INC)

base_target: $(FILES) speedy.o
	$(FC) $(COMPOPTS) $(FILES) speedy.o -o speedy $(LIB)

.PHONY: clean
clean:
	rm -f *.o *.mod

speedy.o               : mod_tsteps.o date.o input_output.o shortwave_radiation.o time_stepping.o\
                         diagnostics.o
boundaries.o           : mod_dyncon1.o mod_atparam.o input_output.o spectral.o
coupler.o              : land_model.o sea_model.o
dyn_geop.o             : mod_atparam.o prognostics.o mod_dyncon1.o
dyn_grtend.o           : mod_atparam.o prognostics.o mod_dyncon1.o mod_dyncon2.o physics.o
dyn_sptend.o           : mod_atparam.o prognostics.o mod_dyncon1.o mod_dyncon2.o
humidity.o             : mod_atparam.o
implicit.o             : mod_atparam.o mod_dyncon0.o mod_dyncon1.o mod_dyncon2.o mod_hdifcon.o\
                         mod_tsteps.o
initialization.o       : coupler.o mod_tsteps.o mod_atparam.o date.o input_output.o time_stepping.o\
                         boundaries.o spectral.o sea_model.o physics.o
ini_indyns.o           : mod_tsteps.o mod_dyncon0.o mod_atparam.o mod_dyncon1.o mod_hdifcon.o\
                         spectral.o
ini_invars.o           : mod_dyncon0.o mod_atparam.o prognostics.o mod_dyncon1.o boundaries.o\
                         diagnostics.o
ini_fordate.o          : mod_dyncon0.o shortwave_radiation.o mod_atparam.o\
 						 mod_hdifcon.o physical_constants.o boundaries.o date.o land_model.o\
						 mod_radcon.o surface_fluxes.o date.o sea_model.o longwave_radiation.o\
						 humidity.o
land_model.o           : mod_atparam.o date.o mod_cpl_bcinterp.o input_output.o\
                         boundaries.o
date.o            	   : mod_tsteps.o
diagnostics.o          : mod_tsteps.o mod_atparam.o prognostics.o spectral.o
mod_dyncon1.o          : mod_atparam.o
mod_dyncon2.o          : mod_atparam.o
prognostics.o          : mod_atparam.o
mod_hdifcon.o          : mod_atparam.o
input_output.o         : mod_atparam.o mod_dyncon1.o physical_constants.o date.o mod_tsteps.o\
                         prognostics.o spectral.o
physical_constants.o   : mod_atparam.o
mod_radcon.o           : mod_atparam.o
physics.o              : mod_atparam.o coupler.o physical_constants.o boundaries.o land_model.o\
                         sea_model.o sppt.o mod_tsteps.o precipitation.o surface_fluxes.o\
						 vertical_diffusion.o shortwave_radiation.o longwave_radiation.o\
						 mod_dyncon1.o humidity.o
longwave_radiation.o   : mod_atparam.o physical_constants.o mod_radcon.o
sea_model.o            : mod_atparam.o input_output.o boundaries.o mod_dyncon1.o mod_cpl_bcinterp.o\
 						 date.o mod_tsteps.o
shortwave_radiation.o  : mod_atparam.o physical_constants.o mod_radcon.o
surface_fluxes.o       : mod_atparam.o physical_constants.o mod_radcon.o land_model.o humidity.o

ppo_dmflux.o           : mod_tsteps.o mod_atparam.o physical_constants.o land_model.o sea_model.o\
                         physics.o mod_radcon.o
precipitation.o        : mod_atparam.o physical_constants.o mod_dyncon1.o
spectral.o             : mod_atparam.o
spe_spectral.o         : mod_atparam.o spectral.o
spe_subfft_fftpack.o   : mod_atparam.o spectral.o
sppt.o                 : mod_atparam.o mod_tsteps.o mod_dyncon1.o spectral.o
time_stepping.o        : mod_tsteps.o mod_dyncon0.o mod_atparam.o prognostics.o\
 					     mod_hdifcon.o implicit.o
vertical_diffusion     : mod_atparam.o physical_constants.o
