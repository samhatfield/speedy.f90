SPEEDY.f90 is an intermediate complexity atmospheric general circulation model written in modern Fortran. It is based on [SPEEDY](http://users.ictp.it/~kucharsk/speedy-net.html), developed by Fred Kucharski, Franco Molteni and Martin P. King.

## Installation

SPEEDY.f90 has only one dependency: the NetCDF library. To build SPEEDY.f90:

1. Install the [NetCDF library](https://www.unidata.ucar.edu/software/netcdf) and locate the `netcdf.mod` file. For example, on my system it is stored in `/usr/include`.
2. Set the `NETCDF` environment variable to point to the directory containing `netcdf.mod`. For example, for my system I run `export NETCDF=/usr/include`.
3. Run `build.sh` to build SPEEDY.f90: `bash build.sh`. A binary directory, `bin`, will be created an the SPEEDY.f90 executable `speedy` will be placed in this directory.
4. Run `run.sh` to run SPEEDY.f90: `bash run.sh`. The output will be stored in `rundir`. By default, SPEEDY.f90 will run for two days and output one NetCDF file for each time step.