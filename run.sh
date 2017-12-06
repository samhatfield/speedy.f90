#!/bin/bash

# $1 = resolution (eg t21, t30)
# $2 = experiment no. (eg 111)
# $3 = experiment no. for restart file ( 0 = no restart ) 
# $4 = make only, and don't run? ("make" for yes, "run" for no)


if [ $# -ne 4 ] ; then
    echo 'Usage: '$0' resol. exp_no. restart_no make_only_or_run' 1>&2
    exit 1
fi

# Define directory names
UT=`pwd`
SRC=$UT/source	
TMP=$UT/tmp
mkdir -p $UT/output/exp_$2	
OUT=$UT/output/exp_$2
CD=$UT/output/exp_$3	

# Copy files from basic version directory

echo "copying from $SRC/source to $TMP"

mkdir -p $TMP
cd $TMP
rm *

cp $SRC/*.f90      $TMP/
cp $SRC/*.h      $TMP/
cp $SRC/*.s      $TMP/
cp $SRC/makefile $TMP/

# Set experiment no. and restart file (if needed)

echo $3 >  fort.2
echo $2 >> fort.2

if [ $3 != 0 ] ; then
  echo "link restart file atgcm$3.rst to fort.3"
  ln -s $CD/atgcm$3.rst fort.3
fi 

# Link input files

echo 'link input files to fortran units'

ksh inpfiles.s $1

ls -l fort.*

echo ' compiling at_gcm - calling make'

make clean
make imp.exe || { echo "Compilation failed"; exit 1; }

if [ $4 == make ] ; then
    exit 0
fi

# Write date input file
cat << EOF > fort.2
0
1982
01
01
00
EOF

time ./imp.exe | tee out.lis

mv out.lis $OUT/atgcm$2.lis
mv fort.10 $OUT/atgcm$2.rst

mv at*$2.ctl   $OUT
mv at*$2_*.grd $OUT

mv day*$2.ctl   $OUT
mv day*$2_*.grd $OUT

mv *.grd $OUT
cp ../ctl_files/yyyymmddhh.ctl $OUT

cd $OUT
