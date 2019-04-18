#!/bin/bash

# Define directory names
ROOT=`pwd`
BIN=$ROOT/bin
RUNDIR=$ROOT/rundir
CLIM=$ROOT/data/bc/t30/clim
ANOM=$ROOT/data/bc/t30/anom

# Check model executable exists
if [ ! -f $BIN/speedy ]; then
    echo "No executable found in bin ($BIN)"
    echo "Have you run ./build.sh yet?"
fi

# Make and move to run directory
rm -rf $RUNDIR
mkdir -p $RUNDIR
cd $RUNDIR

# Copy executable
cp $BIN/speedy .

# Link input files
ln -s $CLIM/sfc.grd   fort.20
ln -s $CLIM/sst.grd   fort.21
ln -s $CLIM/icec.grd  fort.22
ln -s $CLIM/stl.grd   fort.23
ln -s $CLIM/snowd.grd fort.24
ln -s $CLIM/swet.grd  fort.26
ln -s $ANOM/ssta.grd  fort.30

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
time ./speedy | tee output.txt
