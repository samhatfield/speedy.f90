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

# Copy namelist file to run directory
cp ../namelist.nml $RUNDIR

# Run SPEEDY
if [ "$1" = "--profile" ]; then
    if ! [ -x "$(command -v dot)" ]; then
        echo "dot command not found"
        echo "You must install graphviz to use the --profile option"
        exit 1
    fi

    # Run SPEEDY and generate profile data
    ./speedy
    gprof speedy gmon.out > profile.txt

    # Remove __*_MOD_ function prefixes from profiler output
    sed -e "s/ __.*_MOD_//g" -i profile.txt

    # Generate call graph using gprof2dot and graphviz (dot)
    python $ROOT/scripts/gprof2dot.py --skew 0.1 -n 0.8 profile.txt\
        | dot -Tpdf -o profile.pdf
else
    time ./speedy | tee output.txt
fi
