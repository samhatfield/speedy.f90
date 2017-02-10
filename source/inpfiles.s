# Script to link fortran units to input files; $1 = resolution (t21, t30, ..)

SB=../data/bc/$1/clim
SC=../data/bc/$1/anom
SH=../hflux	

 ln -sf $SB/orog_lsm_alb.${1}.grd         fort.20
 ln -sf $SB/sst_7908clim.${1}.sea.grd     fort.21
 ln -sf $SB/seaice_7908clim.${1}.sea.grd  fort.22
 ln -sf $SB/surfv_st3_7908clim.${1}.land.grd   fort.23	
 ln -sf $SB/sndep_7908clim.${1}.land.grd  fort.24
 ln -sf $SB/veget.${1}.land.grd           fort.25
 ln -sf $SB/soilw_7908clim.${1}.land.grd  fort.26

 cp    $SC/noaa_anom_1854_2010_mean1979_2008.${1}.grd fort.30

	
cp    $SH/hflux_speedy_ver41_1979_2008_clim.grd  fort.31
#cp    $SH/hflux_zero.grd  fort.31	

#ln -s $SB/ocean_model_sst.grd  fort.50
