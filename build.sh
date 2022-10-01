#!/bin/sh

# Check that NetCDF variable is set
if [ -z "$NETCDF" ]; then
    echo "\$NETCDF variable not set"
    echo "Aborting build"
    exit 1
fi

#
MAKE=$(uname -a | awk '$1~"BSD"{print "g"}')make

# Name of makefile
MAKEFILE=gfortran.makefile

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
cp $SRC/$MAKEFILE .

# Compile SPEEDY and delete source files
$MAKE -f $MAKEFILE -s clean
echo 'Compiling SPEEDY'
if [ "$1" = "--profile" ]; then
    $MAKE -f $MAKEFILE -s profile || { echo "Compilation failed"; exit 1; }
else
    $MAKE -f $MAKEFILE -s || { echo "Compilation failed"; exit 1; }
fi

rm *.f90 *.o $MAKEFILE *.mod
