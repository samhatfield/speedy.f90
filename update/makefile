###COMPILE=f77
###COMPILE=pgf90 
COMPILE=gfortran


 COMLIB1=

#COMOTT1=-dalign -O3 -r8
###COMOTT1=-dalign -fast -r8
#COMOTT1=-Ofast -mips4 -IPA -r8
# COMOTT1=  -fast -O3 -r8
 COMOTT1=  -Ofast -fdefault-real-8
### COMOTT1=  -fast  

# COMCONV= -byteswapio 
 COMCONV= -fconvert=swap

FILES= \
       cpl_bcinterp.o \
       cpl_land.o \
       cpl_main_interface.o \
       cpl_sea.o \
       cpl_sea_model.o \
       dyn_geop.o \
       dyn_grtend.o \
       dyn_implic.o \
       dyn_sptend.o \
       dyn_step.o \
       dyn_stloop.o \
       ini_agcm_init.o \
       ini_impint.o \
       ini_indyns.o \
       ini_inbcon.o \
       ini_iniatm.o \
       ini_inphys.o \
       ini_inirdf.o \
       ini_invars.o \
       ini_stepone.o \
       ini_fordate.o \
	   mod_atparam.o \
	   mod_cli_land.o \
	   mod_cnvcon.o \
	   mod_cplcon_sea.o \
	   mod_cplvar_sea.o \
	   mod_cpl_flags.o \
       mod_cpl_land_model.o \
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
	   mod_physcon.o \
	   mod_randfor.o \
	   mod_sflcon.o \
	   mod_surfcon.o \
	   mod_tmean.o \
	   mod_tsteps.o \
	   mod_var_land.o \
	   mod_vdicon.o \
       phy_convmf.o \
       phy_lscond.o \
       phy_phypar.o \
       phy_radiat.o \
       phy_shtorh.o \
       phy_suflux.o \
       phy_vdifsc.o \
       ppo_diagns.o \
       ppo_restart.o \
       ppo_setctl.o \
       ppo_setgrd.o \
       ppo_tminc.o \
       ppo_dmflux.o \
       ppo_tmout.o \
       ppo_dmout.o \
       spe_matinv.o \
       spe_spectral.o \
       spe_subfft_fftpack.o

#       land2atm.o \
#       land_init.o \
#       sea2atm.o \
#       sea_init.o \
#       get_from_coupler.o \
#       send_to_coupler.o \

.SUFFIXES: .f .f90 .o

.f90.o:
	$(COMPILE) $(COMCONV) $(COMOTT1) -c $<

.f.o:
	$(COMPILE) $(COMCONV) $(COMOTT1) -c $<

imp.exe: $(FILES)  at_gcm.o
	$(COMPILE) $(COMCONV) $(COMOTT1) $(FILES) at_gcm.o -o imp.exe $(COMLIB1) 


at_gcm.o               : at_gcm.f90 mod_tsteps.o mod_date.o
#addflx.o               : addflx.f
#atm2land.o             : atm2land.f
#atm2sea.o              : atm2sea.f
#get_from_coupler.o     : get_from_coupler.f
cpl_bcinterp.o         : cpl_bcinterp.f90
cpl_land.o             : cpl_land.f90 mod_cpl_flags.o mod_atparam.o mod_cpl_land_model.o mod_flx_land.o\
						 mod_cli_land.o mod_date.o mod_var_land.o
mod_cpl_land_model.o   : mod_cpl_land_model.f90 mod_atparam.o
cpl_main_interface.o   : cpl_main_interface.f90 mod_atparam.o mod_cpl_land_model.o mod_surfcon.o mod_cli_land.o
cpl_sea.o              : cpl_sea.f90 mod_cpl_flags.o mod_atparam.o mod_cplvar_sea.o mod_date.o mod_flx_sea.o
cpl_sea_model.o        : cpl_sea_model.f90 mod_atparam.o mod_cplcon_sea.o mod_cplvar_sea.o
dyn_geop.o             : dyn_geop.f90 mod_atparam.o mod_dynvar.o mod_dyncon1.o
dyn_grtend.o           : dyn_grtend.f90 mod_atparam.o mod_dynvar.o mod_dyncon1.o mod_dyncon2.o
dyn_implic.o           : dyn_implic.f90 mod_atparam.o mod_dyncon1.o mod_dyncon1.o mod_dyncon2.o
dyn_sptend.o           : dyn_sptend.f90 mod_atparam.o mod_dynvar.o mod_dyncon1.o mod_dyncon2.o
dyn_step.o             : dyn_step.f90 mod_tsteps.o mod_dyncon0.o mod_atparam.o mod_dynvar.o mod_hdifcon.o
dyn_stloop.o           : dyn_stloop.f90 mod_lflags.o mod_tsteps.o
#ini_coupler.o          : ini_coupler.f
ini_agcm_init.o        : ini_agcm_init.f90 mod_cpl_flags.o mod_tsteps.o mod_atparam.o mod_date.o
ini_impint.o           : ini_impint.f90 mod_dyncon0.o mod_atparam.o mod_dyncon1.o mod_dyncon2.o mod_hdifcon.o
ini_indyns.o           : ini_indyns.f90 mod_tsteps.o mod_dyncon0.o mod_atparam.o mod_dyncon1.o mod_hdifcon.o
ini_inbcon.o           : ini_inbcon.f90 mod_cpl_flags.o mod_tsteps.o mod_dyncon0.o mod_atparam.o mod_dyncon1.o\
						 mod_surfcon.o mod_cli_land.o
ini_iniatm.o           : ini_iniatm.f90 mod_tsteps.o mod_atparam.o mod_dyncon1.o mod_tmean.o mod_date.o
ini_inphys.o           : ini_inphys.f90 mod_atparam.o mod_physcon.o
ini_inirdf.o           : ini_inirdf.f90 mod_atparam.o mod_physcon.o mod_randfor.o
ini_invars.o           : ini_invars.f90 mod_tsteps.o mod_dyncon0.o mod_atparam.o mod_dynvar.o mod_dyncon1.o\
						 mod_surfcon.o mod_date.o
ini_stepone.o          : ini_stepone.f90 mod_tsteps.o
ini_fordate.o          : ini_fordate.f90 mod_dyncon0.o mod_lflags.o mod_atparam.o mod_hdifcon.o mod_physcon.o\
						 mod_surfcon.o mod_cli_land.o mod_date.o mod_var_land.o
#land2atm.o             : land2atm.f
#land_init.o            : land_init.f
mod_atparam.o          : mod_atparam.f90
mod_cli_land.o         : mod_cli_land.f90 mod_atparam.o
mod_cnvcon.o           : mod_cnvcon.f90
mod_cplcon_sea.o       : mod_cplcon_sea.f90 mod_atparam.o
mod_cplvar_sea.o       : mod_cplvar_sea.f90 mod_atparam.o
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
mod_physcon.o          : mod_physcon.f90 mod_atparam.o
mod_randfor.o          : mod_randfor.f90 mod_atparam.o
mod_sflcon.o           : mod_sflcon.f90 mod_atparam.o
mod_surfcon.o          : mod_surfcon.f90 mod_atparam.o
mod_tmean.o            : mod_tmean.f90 mod_atparam.o
mod_tsteps.o           : mod_tsteps.f90
mod_var_land.o         : mod_var_land.f90 mod_atparam.o
mod_vdicon.o           : mod_vdicon.f90
phy_convmf.o           : phy_convmf.f90 mod_cnvcon.o mod_atparam.o mod_physcon.o
phy_lscond.o           : phy_lscond.f90 mod_lsccon.o mod_atparam.o mod_physcon.o
phy_phypar.o           : phy_phypar.f90 mod_cpl_flags.o mod_lflags.o mod_atparam.o mod_physcon.o mod_randfor.o\
					     mod_surfcon.o mod_var_land.o
phy_radiat.o           : phy_radiat.f90 mod_atparam.o mod_physcon.o
phy_shtorh.o           : phy_shtorh.f90
phy_suflux.o           : phy_suflux.f90 mod_atparam.o mod_sflcon.o mod_physcon.o
phy_vdifsc.o           : phy_vdifsc.f90 mod_atparam.o mod_vdicon.o mod_physcon.o
ppo_diagns.o           : ppo_diagns.f90 mod_tsteps.o mod_atparam.o mod_dynvar.o
ppo_restart.o          : ppo_restart.f90 mod_tsteps.o mod_atparam.o mod_dynvar.o mod_date.o
ppo_setctl.o           : ppo_setctl.f90
ppo_setgrd.o           : ppo_setgrd.f90 mod_tsteps.o mod_date.o
ppo_tminc.o            : ppo_tminc.f90 mod_lflags.o mod_atparam.o mod_tmean.o mod_physcon.o mod_surfcon.o\
						 mod_cli_land.o
ppo_dmflux.o           : ppo_dmflux.f90 mod_tsteps.o mod_atparam.o mod_tmean.o mod_flx_land.o mod_physcon.o\
						 mod_surfcon.o mod_flx_sea.o
ppo_tmout.o            : ppo_tmout.f90 mod_tsteps.o mod_atparam.o mod_tmean.o mod_physcon.o
ppo_dmout.o            : ppo_dmout.f90 mod_tsteps.o mod_atparam.o mod_tmean.o
#sea2atm.o              : sea2atm.f
#sea_init.o             : sea_init.f 
#send_to_coupler.o      : send_to_coupler.f
spe_matinv.o           : spe_matinv.f90
spe_spectral.o         : spe_spectral.f90 mod_atparam.o
spe_subfft_fftpack.o   : spe_subfft_fftpack.f90 mod_atparam.o
