# Fortran compiler
FC=gfortran

# Set variables depending on target
speedy : COMPOPTS = $(OPT) $(BASE)
profile : COMPOPTS = -pg $(BASE)

# Default
speedy : main

# For profiling
profile : main

# Base compiler options (always used)
BASE=-fdefault-real-8 -fconvert=swap

# Optimisation flags (disabled for debugging, profiling etc.)
OPT=-Ofast

# # Location of NetCDF module (netcdf.mod)
INC=-I$(NETCDF)

# Library flags
LIB=-lnetcdff -lnetcdf

FILES= \
       cpl_main_interface.o \
       dyn_geop.o \
       dyn_grtend.o \
       dyn_implic.o \
       dyn_sptend.o \
       dyn_step.o \
       ini_agcm_init.o \
       ini_impint.o \
       ini_indyns.o \
       ini_inbcon.o \
       ini_iniatm.o \
       ini_inphys.o \
       ini_invars.o \
       ini_stepone.o \
       ini_fordate.o \
	   mod_atparam.o \
	   mod_cli_sea.o \
	   mod_cnvcon.o \
	   mod_cpl_bcinterp.o \
	   mod_cpl_flags.o \
       mod_cpl_land_model.o \
       mod_cpl_sea_model.o \
	   mod_date.o \
	   mod_dyncon0.o \
	   mod_dyncon1.o \
	   mod_dyncon2.o \
	   mod_dynvar.o \
	   mod_flx_land.o \
	   mod_flx_sea.o \
	   mod_hdifcon.o \
	   mod_lflags.o \
	   mod_lsccon.o \
	   mod_input.o \
	   mod_output.o \
	   mod_physcon.o \
	   mod_physvar.o \
	   mod_radcon.o \
	   mod_sflcon.o \
	   mod_spectral.o \
	   mod_sppt.o \
	   mod_surfcon.o \
	   mod_tsteps.o \
	   mod_var_sea.o \
	   mod_vdicon.o \
       phy_convmf.o \
       phy_lscond.o \
       phy_phypar.o \
       phy_radiat.o \
       phy_shtorh.o \
       phy_suflux.o \
       phy_vdifsc.o \
       ppo_diagns.o \
       ppo_dmflux.o \
       spe_matinv.o \
       spe_spectral.o \
       spe_subfft_fftpack.o

%.o: %.f90
	$(FC) $(COMPOPTS) -c $< $(INC)

main: $(FILES) at_gcm.o
	$(FC) $(COMPOPTS) $(FILES) at_gcm.o -o speedy $(LIB)

.PHONY: clean
clean:
	rm -f *.o *.mod

at_gcm.o               : at_gcm.f90 mod_tsteps.o mod_date.o
mod_cpl_land_model.o   : mod_cpl_land_model.f90 mod_atparam.o mod_cpl_flags.o mod_flx_land.o\
                         mod_date.o mod_cpl_bcinterp.o
cpl_main_interface.o   : cpl_main_interface.f90 mod_atparam.o mod_cpl_land_model.o\
 						 mod_cpl_sea_model.o mod_surfcon.o mod_cli_sea.o
cpl_sea.o              : cpl_sea.f90 mod_cpl_flags.o mod_atparam.o mod_date.o\
 						 mod_flx_sea.o mod_cli_sea.o mod_var_sea.o mod_input.o
mod_cpl_sea_model.o    : mod_cpl_sea_model.f90 mod_atparam.o\
                         mod_cli_sea.o mod_var_sea.o mod_flx_sea.o mod_input.o
dyn_geop.o             : dyn_geop.f90 mod_atparam.o mod_dynvar.o mod_dyncon1.o
dyn_grtend.o           : dyn_grtend.f90 mod_atparam.o mod_dynvar.o mod_dyncon1.o mod_dyncon2.o
dyn_implic.o           : dyn_implic.f90 mod_atparam.o mod_dyncon1.o mod_dyncon1.o mod_dyncon2.o
dyn_sptend.o           : dyn_sptend.f90 mod_atparam.o mod_dynvar.o mod_dyncon1.o mod_dyncon2.o
dyn_step.o             : dyn_step.f90 mod_tsteps.o mod_dyncon0.o mod_atparam.o mod_dynvar.o\
 					     mod_hdifcon.o
ini_agcm_init.o        : ini_agcm_init.f90 mod_cpl_flags.o mod_tsteps.o mod_atparam.o mod_date.o
ini_impint.o           : ini_impint.f90 mod_dyncon0.o mod_atparam.o mod_dyncon1.o mod_dyncon2.o\
 					     mod_hdifcon.o
ini_indyns.o           : ini_indyns.f90 mod_tsteps.o mod_dyncon0.o mod_atparam.o mod_dyncon1.o\
 						 mod_hdifcon.o mod_spectral.o
ini_inbcon.o           : ini_inbcon.f90 mod_cpl_flags.o mod_tsteps.o mod_dyncon0.o mod_atparam.o\
                         mod_dyncon1.o mod_surfcon.o mod_cli_sea.o mod_input.o
ini_iniatm.o           : ini_iniatm.f90 mod_tsteps.o mod_atparam.o mod_dyncon1.o mod_date.o\
 						 mod_spectral.o mod_output.o
ini_inphys.o           : ini_inphys.f90 mod_atparam.o mod_physcon.o
ini_invars.o           : ini_invars.f90 mod_tsteps.o mod_dyncon0.o mod_atparam.o mod_dynvar.o\
                         mod_dyncon1.o mod_surfcon.o mod_date.o
ini_stepone.o          : ini_stepone.f90 mod_tsteps.o
ini_fordate.o          : ini_fordate.f90 mod_dyncon0.o mod_lflags.o mod_atparam.o mod_hdifcon.o\
                         mod_physcon.o mod_surfcon.o mod_date.o mod_cpl_land_model.o\
						 mod_cli_sea.o mod_var_sea.o mod_radcon.o
mod_atparam.o          : mod_atparam.f90
mod_cli_sea.o          : mod_cli_sea.f90 mod_atparam.o
mod_cnvcon.o           : mod_cnvcon.f90
mod_cpl_bcinterp.o     : mod_cpl_bcinterp.f90
mod_cpl_flags.o        : mod_cpl_flags.f90
mod_date.o             : mod_date.f90 mod_atparam.o mod_tsteps.o
mod_dyncon0.o          : mod_dyncon0.f90
mod_dyncon1.o          : mod_dyncon1.f90 mod_atparam.o
mod_dyncon2.o          : mod_dyncon2.f90 mod_atparam.o
mod_dynvar.o           : mod_dynvar.f90 mod_atparam.o
mod_flx_land.o         : mod_flx_land.f90 mod_atparam.o
mod_flx_sea.o          : mod_flx_sea.f90 mod_atparam.o
mod_hdifcon.o          : mod_hdifcon.f90 mod_atparam.o
mod_lflags.o           : mod_lflags.f90
mod_lsccon.o           : mod_lsccon.f90
mod_input.o            : mod_input.f90 mod_atparam.o
mod_output.o           : mod_output.f90 mod_atparam.o mod_dyncon1.o mod_physcon.o mod_date.o\
                         mod_tsteps.o mod_dynvar.o
mod_physcon.o          : mod_physcon.f90 mod_atparam.o
mod_physvar.o          : mod_physvar.f90 mod_atparam.o
mod_radcon.o           : mod_radcon.f90 mod_atparam.o
mod_sflcon.o           : mod_sflcon.f90 mod_atparam.o
mod_spectral.o         : mod_spectral.f90 mod_atparam.o
mod_sppt.o             : mod_atparam.o mod_tsteps.o mod_dyncon1.o mod_spectral.o
mod_surfcon.o          : mod_surfcon.f90 mod_atparam.o
mod_tsteps.o           : mod_tsteps.f90
mod_var_sea.o          : mod_var_sea.f90 mod_atparam.o
mod_vdicon.o           : mod_vdicon.f90
phy_convmf.o           : phy_convmf.f90 mod_cnvcon.o mod_atparam.o mod_physcon.o
phy_lscond.o           : phy_lscond.f90 mod_lsccon.o mod_atparam.o mod_physcon.o
phy_phypar.o           : phy_phypar.f90 mod_cpl_flags.o mod_lflags.o mod_atparam.o mod_physcon.o\
					     mod_surfcon.o mod_cpl_land_model.o mod_cli_sea.o mod_var_sea.o\
						 mod_physvar.o mod_sppt.o mod_tsteps.o
phy_radiat.o           : phy_radiat.f90 mod_atparam.o mod_physcon.o mod_radcon.o
phy_shtorh.o           : phy_shtorh.f90
phy_suflux.o           : phy_suflux.f90 mod_atparam.o mod_sflcon.o mod_physcon.o mod_radcon.o
phy_vdifsc.o           : phy_vdifsc.f90 mod_atparam.o mod_vdicon.o mod_physcon.o
ppo_diagns.o           : ppo_diagns.f90 mod_tsteps.o mod_atparam.o mod_dynvar.o mod_date.o
ppo_dmflux.o           : ppo_dmflux.f90 mod_tsteps.o mod_atparam.o mod_flx_land.o mod_physcon.o\
						 mod_surfcon.o mod_flx_sea.o mod_var_sea.o mod_physvar.o mod_radcon.o
spe_matinv.o           : spe_matinv.f90
spe_spectral.o         : spe_spectral.f90 mod_atparam.o mod_spectral.o
spe_subfft_fftpack.o   : spe_subfft_fftpack.f90 mod_atparam.o mod_spectral.o
