###COMPILE=f77
COMPILE=pgf90 


 COMLIB1=

#COMOTT1=-dalign -O3 -r8
###COMOTT1=-dalign -fast -r8
#COMOTT1=-Ofast -mips4 -IPA -r8
 COMOTT1=  -fast -O3 -r8
### COMOTT1=  -fast  

 COMCONV= -byteswapio 

FILES= \
       cpl_bcinterp.o \
       cpl_land.o \
       cpl_land_model.o \
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

.SUFFIXES: .f .o

.f.o:
	$(COMPILE) $(COMCONV) $(COMOTT1) -c $<

imp.exe: $(FILES)  at_gcm.o
	$(COMPILE) $(COMCONV) $(COMOTT1) $(FILES) at_gcm.o -o imp.exe $(COMLIB1) 


at_gcm.o               : at_gcm.f
#addflx.o               : addflx.f
#atm2land.o             : atm2land.f
#atm2sea.o              : atm2sea.f
#get_from_coupler.o     : get_from_coupler.f
cpl_bcinterp.o         : cpl_bcinterp.f
cpl_land.o             : cpl_land.f
cpl_land_model.o       : cpl_land_model.f
cpl_main_interface.o   : cpl_main_interface.f
cpl_sea.o              : cpl_sea.f
cpl_sea_model.o        : cpl_sea_model.f
dyn_geop.o             : dyn_geop.f
dyn_grtend.o           : dyn_grtend.f
dyn_implic.o           : dyn_implic.f
dyn_sptend.o           : dyn_sptend.f
dyn_step.o             : dyn_step.f
dyn_stloop.o           : dyn_stloop.f
#ini_coupler.o          : ini_coupler.f
ini_agcm_init.o        : ini_agcm_init.f
ini_impint.o           : ini_impint.f
ini_indyns.o           : ini_indyns.f
ini_inbcon.o           : ini_inbcon.f
ini_iniatm.o           : ini_iniatm.f
ini_inphys.o           : ini_inphys.f
ini_inirdf.o           : ini_inirdf.f
ini_invars.o           : ini_invars.f
ini_stepone.o          : ini_stepone.f
ini_fordate.o          : ini_fordate.f
#land2atm.o             : land2atm.f
#land_init.o            : land_init.f
phy_convmf.o           : phy_convmf.f
phy_lscond.o           : phy_lscond.f
phy_phypar.o           : phy_phypar.f
phy_radiat.o           : phy_radiat.f
phy_shtorh.o           : phy_shtorh.f
phy_suflux.o           : phy_suflux.f
phy_vdifsc.o           : phy_vdifsc.f
ppo_diagns.o           : ppo_diagns.f
ppo_restart.o          : ppo_restart.f
ppo_setctl.o           : ppo_setctl.f
ppo_setgrd.o           : ppo_setgrd.f
ppo_tminc.o            : ppo_tminc.f
ppo_dmflux.o           : ppo_dmflux.f
ppo_tmout.o            : ppo_tmout.f
ppo_dmout.o            : ppo_dmout.f
#sea2atm.o              : sea2atm.f
#sea_init.o             : sea_init.f 
#send_to_coupler.o      : send_to_coupler.f
spe_matinv.o           : spe_matinv.f
spe_spectral.o         : spe_spectral.f
spe_subfft_fftpack.o   : spe_subfft_fftpack.f 
