![Zenodo DOI: 10.5281/zenodo.5816982](https://zenodo.org/badge/DOI/10.5281/zenodo.5816982.svg)

speedy.f90 is an intermediate complexity atmospheric general circulation model written in modern Fortran. It is based on [SPEEDY](http://users.ictp.it/~kucharsk/speedy-net.html), developed by Fred Kucharski, Franco Molteni and Martin P. King.

## Installation

speedy.f90 has only one dependency: the NetCDF library. To build speedy.f90:

1. Install the [NetCDF library](https://www.unidata.ucar.edu/software/netcdf) and locate the install directory. For example, on my system it is stored in `/usr`.
2. Set the `NETCDF` environment variable to point to the directory containing the NetCDF `include` and `lib` directories. For example, for my system I run `export NETCDF=/usr`.
3. Run `build.sh` to build speedy.f90: `bash build.sh`. A binary directory, `bin`, will be created an the speedy.f90 executable `speedy` will be placed in this directory.
4. Run `run.sh` to run speedy.f90: `bash run.sh`. The output will be stored in `rundir`. By default, speedy.f90 will run for two days and output one NetCDF file for each time step.

## How to Cite

Cite the Zenodo DOI given above directly, as well as the following three model description and verification papers based on the original SPEEDY model:

- Molteni, F., Atmospheric simulations using a GCM with simplified physical parametrizations. I: model climatology and variability in multi-decadal experiments. Climate Dynamics 20, 175–191 (2003). https://doi.org/10.1007/s00382-002-0268-2
- Kucharski, F., Molteni, F. & Bracco, A. Decadal interactions between the western tropical Pacific and the North Atlantic Oscillation. Climate Dynamics 26, 79–91 (2006). https://doi.org/10.1007/s00382-005-0085-5
- Kucharski, F., Molteni, F., King, M. P., Farneti, R., Kang, I., & Feudale, L., On the Need of Intermediate Complexity General Circulation Models: A “SPEEDY” Example. Bulletin of the American Meteorological Society, 94(1), 25-30 (2013). https://doi.org/10.1175/BAMS-D-11-00238.1
