#!/bin/bash

# Check that NetCDF variable is set
if [ -z "$NETCDF" ]; then
    echo "\$NETCDF variable not set"
    echo "Aborting build"
    exit 1
fi

# Define directory names
ROOT=`pwd`
SRC=$ROOT/source
BIN=$ROOT/bin

# Make binary directory
rm -rf $BIN
mkdir $BIN
cd $BIN

# Copy source files
cp $SRC/*.f90    .
cp $SRC/makefile .

# Compile SPEEDY and delete source files
make -s clean
echo 'Compiling SPEEDY'
if [ "$1" = "--profile" ]; then
    make -s profile || { echo "Compilation failed"; exit 1; }
else
    make -s || { echo "Compilation failed"; exit 1; }
fi

rm *.f90 *.o makefile *.mod
