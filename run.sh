#!/bin/bash

# Define directory names
ROOT=`pwd`
SRC=$ROOT/source
RUNDIR=$ROOT/rundir
CLIM=$ROOT/data/bc/t30/clim
ANOM=$ROOT/data/bc/t30/anom

# Make run directory
rm -rf $RUNDIR
mkdir -p $RUNDIR
cd $RUNDIR

# Copy source files
cp $SRC/*.f90    $RUNDIR
cp $SRC/*.h      $RUNDIR
cp $SRC/makefile $RUNDIR

# Link input files
ln -s $CLIM/sfc.grd   fort.20
ln -s $CLIM/sst.grd   fort.21
ln -s $CLIM/icec.grd  fort.22
ln -s $CLIM/stl.grd   fort.23
ln -s $CLIM/snowd.grd fort.24
ln -s $CLIM/swet.grd  fort.26
ln -s $ANOM/ssta.grd  fort.30

# Compile SPEEDY and delete source files
echo 'Compiling SPEEDY'
make -s clean
make -s imp.exe || { echo "Compilation failed"; exit 1; }
rm *.f90 *.h *.o makefile *.mod

# Write date input file
cat << EOF >> fort.2
1982
01
01
00
00
1982
01
02
00
00
EOF

# Run SPEEDY
time ./imp.exe | tee output.txt
