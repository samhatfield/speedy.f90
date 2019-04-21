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
       coupler.o \
	   date.o \
       dyn_geop.o \
       dyn_grtend.o \
       dyn_implic.o \
       dyn_sptend.o \
       dyn_step.o \
	   humidity.o \
       ini_agcm_init.o \
       ini_impint.o \
       ini_indyns.o \
       ini_inbcon.o \
       ini_iniatm.o \
       ini_invars.o \
       ini_stepone.o \
       ini_fordate.o \
	   land_model.o \
	   mod_atparam.o \
	   mod_cpl_bcinterp.o \
	   mod_dyncon0.o \
	   mod_dyncon1.o \
	   mod_dyncon2.o \
	   mod_hdifcon.o \
	   mod_input.o \
	   mod_output.o \
	   mod_physcon.o \
	   mod_radcon.o \
	   mod_spectral.o \
	   mod_surfcon.o \
	   mod_tsteps.o \
	   physics.o \
       phy_radiat.o \
       ppo_diagns.o \
       ppo_dmflux.o \
	   precipitation.o \
	   prognostics.o \
	   sea_model.o \
	   shortwave_radiation.o \
       spe_matinv.o \
       spe_spectral.o \
       spe_subfft_fftpack.o \
	   sppt.o \
	   surface_fluxes.o \
	   vertical_diffusion.o

%.o: %.f90
	$(FC) $(COMPOPTS) -c $< $(INC)

base_target: $(FILES) speedy.o
	$(FC) $(COMPOPTS) $(FILES) speedy.o -o speedy $(LIB)

.PHONY: clean
clean:
	rm -f *.o *.mod

speedy.o               : speedy.f90 mod_tsteps.o date.o mod_output.o shortwave_radiation.o
coupler.o              : coupler.f90 mod_atparam.o land_model.o sea_model.o mod_surfcon.o
cpl_sea.o              : cpl_sea.f90 coupler.o mod_atparam.o date.o mod_input.o
dyn_geop.o             : dyn_geop.f90 mod_atparam.o prognostics.o mod_dyncon1.o
dyn_grtend.o           : dyn_grtend.f90 mod_atparam.o prognostics.o mod_dyncon1.o mod_dyncon2.o\
                         physics.o
dyn_implic.o           : dyn_implic.f90 mod_atparam.o mod_dyncon1.o mod_dyncon1.o mod_dyncon2.o
dyn_sptend.o           : dyn_sptend.f90 mod_atparam.o prognostics.o mod_dyncon1.o mod_dyncon2.o
dyn_step.o             : dyn_step.f90 mod_tsteps.o mod_dyncon0.o mod_atparam.o prognostics.o\
 					     mod_hdifcon.o
humidity.o             : humidity.f90 mod_atparam.o
ini_agcm_init.o        : ini_agcm_init.f90 coupler.o mod_tsteps.o mod_atparam.o date.o
ini_impint.o           : ini_impint.f90 mod_dyncon0.o mod_atparam.o mod_dyncon1.o mod_dyncon2.o\
 					     mod_hdifcon.o
ini_indyns.o           : ini_indyns.f90 mod_tsteps.o mod_dyncon0.o mod_atparam.o mod_dyncon1.o\
 						 mod_hdifcon.o mod_spectral.o
ini_inbcon.o           : ini_inbcon.f90 coupler.o mod_tsteps.o mod_dyncon0.o mod_atparam.o\
                         mod_dyncon1.o mod_surfcon.o mod_input.o
ini_iniatm.o           : ini_iniatm.f90 mod_tsteps.o mod_atparam.o mod_dyncon1.o date.o\
 						 mod_spectral.o mod_output.o
ini_invars.o           : ini_invars.f90 mod_tsteps.o mod_dyncon0.o mod_atparam.o prognostics.o\
                         mod_dyncon1.o mod_surfcon.o date.o
ini_stepone.o          : ini_stepone.f90 mod_tsteps.o
ini_fordate.o          : ini_fordate.f90 mod_dyncon0.o shortwave_radiation.o mod_atparam.o\
 						 mod_hdifcon.o mod_physcon.o mod_surfcon.o date.o land_model.o\
						 mod_radcon.o surface_fluxes.o
land_model.o           : land_model.f90 mod_atparam.o date.o mod_cpl_bcinterp.o
mod_atparam.o          : mod_atparam.f90
mod_cpl_bcinterp.o     : mod_cpl_bcinterp.f90
date.o            	   : date.f90 mod_atparam.o mod_tsteps.o
mod_dyncon0.o          : mod_dyncon0.f90
mod_dyncon1.o          : mod_dyncon1.f90 mod_atparam.o
mod_dyncon2.o          : mod_dyncon2.f90 mod_atparam.o
prognostics.o          : prognostics.f90 mod_atparam.o
mod_hdifcon.o          : mod_hdifcon.f90 mod_atparam.o
mod_input.o            : mod_input.f90 mod_atparam.o
mod_output.o           : mod_output.f90 mod_atparam.o mod_dyncon1.o mod_physcon.o date.o\
                         mod_tsteps.o prognostics.o
mod_physcon.o          : mod_physcon.f90 mod_atparam.o
mod_radcon.o           : mod_radcon.f90 mod_atparam.o
mod_spectral.o         : mod_spectral.f90 mod_atparam.o
mod_surfcon.o          : mod_surfcon.f90 mod_atparam.o
mod_tsteps.o           : mod_tsteps.f90
physics.o              : physics.f90 mod_atparam.o coupler.o mod_physcon.o\
					     mod_surfcon.o land_model.o sea_model.o sppt.o mod_tsteps.o\
						 precipitation.o surface_fluxes.o vertical_diffusion.o shortwave_radiation.o
sea_model.o            : sea_model.f90 mod_atparam.o mod_input.o
phy_radiat.o           : phy_radiat.f90 mod_atparam.o mod_physcon.o mod_radcon.o
shortwave_radiation.o  : shortwave_radiation.f90 mod_atparam.o mod_physcon.o mod_radcon.o
surface_fluxes.o       : surface_fluxes.f90 mod_atparam.o mod_physcon.o mod_radcon.o humidity.o
ppo_diagns.o           : ppo_diagns.f90 mod_tsteps.o mod_atparam.o prognostics.o date.o
ppo_dmflux.o           : ppo_dmflux.f90 mod_tsteps.o mod_atparam.o mod_physcon.o mod_surfcon.o\
                         mod_radcon.o
precipitation.o        : precipitation.f90 mod_atparam.o mod_physcon.o
spe_matinv.o           : spe_matinv.f90
spe_spectral.o         : spe_spectral.f90 mod_atparam.o mod_spectral.o
spe_subfft_fftpack.o   : spe_subfft_fftpack.f90 mod_atparam.o mod_spectral.o
sppt.o                 : sppt.f90 mod_atparam.o mod_tsteps.o mod_dyncon1.o mod_spectral.o
vertical_diffusion     : phy_vdifsc.f90 mod_atparam.o mod_physcon.o
