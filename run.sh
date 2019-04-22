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
ln -s $CLIM/surface.nc .
ln -s $CLIM/sea_surface_temperature.nc .
ln -s $CLIM/sea_ice.nc .
ln -s $CLIM/land.nc .
ln -s $CLIM/snow.nc .
ln -s $CLIM/soil.nc .
ln -s $ANOM/sea_surface_temperature_anomaly.nc .

# Write date input file
cat << EOF >> fort.2
1982
01
01
00
00
1982
01
10
00
00
EOF

# Run SPEEDY
time ./speedy | tee output.txt
