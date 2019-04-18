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
cp $SRC/*.h      .
cp $SRC/makefile .

# Compile SPEEDY and delete source files
echo 'Compiling SPEEDY'
make -s clean
make -s speedy || { echo "Compilation failed"; exit 1; }
rm *.f90 *.h *.o makefile *.mod
