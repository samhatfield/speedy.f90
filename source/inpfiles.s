# Script to link fortran units to input files; $1 = resolution (t21, t30, ..)

SB=../data/bc/$1/clim
SC=../data/bc/$1/anom
SH=../hflux	

ln -sf $SB/sfc.grd   fort.20
ln -sf $SB/sst.grd   fort.21
ln -sf $SB/icec.grd  fort.22
ln -sf $SB/stl.grd   fort.23	
ln -sf $SB/snowd.grd fort.24
ln -sf $SB/swet.grd  fort.26

cp    $SC/ssta.grd  fort.30	
