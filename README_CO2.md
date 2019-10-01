<<<<<<< HEAD
# CO2 Only GCHP on Pleaides
## Prereqs

Pleiades already has almost all the necessary packages available in its module system. To load the appropriate modules, place the following text in a file called  `gchp.ifort18_sgimpi_pleiades.env` in your home directory:
```
#!/bin/bash
# Based on seastham's home/pleiades.basrch

# load ifort, sgi mpi, and netcdf
module load comp-intel/2018.3.222
module load mpi-sgi/mpt
module load hdf4/4.2.12 # required by netcdf
module load hdf5/1.8.18_mpt # also required by netcdf?
module load netcdf/4.4.1.1_mpt


export ESMF_COMM=mpi
export ESMF_COMPILER=intel

# Tell GCHP where the MPI binaries, libraries and so on can be found
export MPI_ROOT=$( dirname $( dirname $( which mpiexec ) ) )

# Set up the compilers
export FC=ifort
export F90=$FC
export F9X=$FC
export F77=$FC
export CC=gcc
export CXX=g++

export OMP_STACKSIZE=500m
ulimit -s unlimited

# Needed one NetCDF is installed
export NETCDF_HOME=$(nc-config --prefix)

export GC_BIN="$NETCDF_HOME/bin"
export GC_INCLUDE="$NETCDF_HOME/include"
export GC_LIB="$NETCDF_HOME/lib"

# If using NetCDF after the C/Fortran split (4.3+), then you will need to
# specify the following additional environment variables
export NETCDF_FORTRAN_HOME=$(nf-config --prefix)
export GC_F_BIN="$NETCDF_FORTRAN_HOME/bin"
export GC_F_INCLUDE="$NETCDF_FORTRAN_HOME/include"
export GC_F_LIB="$NETCDF_FORTRAN_HOME/lib"

export PATH=${NETCDF_FORTRAN_HOME}/bin:$PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${NETCDF_FORTRAN_HOME}/lib

# Add NetCDF to path
export PATH=$PATH:${NETCDF_HOME}/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${NETCDF_HOME}/lib

# Disable OpenMP
export OMP_NUM_THREADS=1

# Attempt an alias for mpifort if not found
#command -v mpifort >/dev/null 2>&1 || alias mpifort=mpif90
export PATH=${PATH}:${HOME}/mpi_extra

# Set path to GMAO Fortran template library (gFTL)
export gFTL=$(readlink -f ./gFTL)

# Set ESMF optimization (g=debugging, O=optimized (capital o))
export ESMF_BOPT=O
```

On library that is not available is the GMOA Fortran library. If you already have > v1.0 available, make sure you know the path it was installed at as you will need it. If you do not, follow these instructions to install it, and, again, make note of the location where you installed it.

### Installing gFTL

1. Navigate to directory where you want to download gFTL
2. Type the following at the command prompt:
   ```
   $ git clone https://github.com/Goddard-Fortran-Ecosystem/gFTL`
   $ cd gFTL
   $ git checkout v1.0.0
   $ cmake . -DCMAKE_INSTALL_PREFIX=.
   $ make install
   ```
3. Verify success by checking that include/templates and include/types exist

## Getting the code
First, checkout the dev/v12.5.0 branch of both the https://github.com/TerribleNews/geos-chem, and inside that folder https://github.com/TerribleNews/gchp.

```
$ git clone https://github.com/TerribleNews/geos-chem.git Code.GCHP_CO2
$ cd Code.GCHP_CO2
$ git checkout --track origin/dev/12.5.0
$ git clone https://github.com/TerribleNews/gchp.git GCHP
$ cd GCHP
$ git checkout --track origin/dev/12.5.0
```


## Building the code
GEOS-Chem and GCHP are built from a run directory. To create a run directory, use the `createRunDirectory.sh` script in the GCHP/Run directory. In the following block, anything you are supposed to enter is prefixed by a `$`
```
$ cd Run
Define paths to ExtData and the Goddard Fortran Template Library (gFTL).
These will be stored in $HOME/.geoschem/config for future automatic use.

Enter path for ExtData:
```
Enter `/nobackupp2/seastham/ExtData` as the path for ExtData

```
$ /nobackupp2/seastham/ExtData

If you have not downloaded gFTL then enter q to exit.
Follow these instructions at the command prompt to install:

      1. Navigate to directory where you want to download gFTL
      2. Type the following at the command prompt:
         $ git clone https://github.com/Goddard-Fortran-Ecosystem/gFTL
         $ cd gFTL
         $ git checkout v1.0.0
         $ cmake . -DCMAKE_INSTALL_PREFIX=.
         $ make install
      3. Verify success by checking that include/templates and include/types exist

Enter path for gFTL:
```
Here you'll need the path to your gFTL install as mentioned before
```
$ /home1/clee59/src/gFTL

Choose simulation type:
  1. TransportTracers
  2. Standard
  3. Benchmark
  4. CO2
$ 4

Choose meteorology source:
  1. GEOS-FP
  2. MERRA2
```
Only the MERRA2 files are currently in ExtData.
```
$ 2

Enter path where the run directory will be created:
```
Feel free to choose your own location to save your runs. It should probably be on the /nobackup drives.
```
$ /nobackupp13/clee59/GCHP_runs/

Enter run directory name, or press return to use default:
$ gchp_co2

Do you want to track run directory changes with git? (y/n)
n
```

## Running GCHP
Go to the directory you created, gchp_co2 and copy the correct restart file and CO2 emissions:
```
$ rm initial_GEOSChem_rst.c*.nc
$ cp /nobackupp13/clee59/SHARED/* ./
```

Create a link to the environment file you created at the beginning of this README:
```
$ ./setEnvironment ~/gchp.ifort18_sgimpi_pleiades.env
```

Then copy and edit the run script:
```
$ cp runScriptSamples/gchp.pleiades.run ./
```

Inside `gchp.pleiades.run` change the email address and job account group at the top of the script to match your email and job account
```
#PBS -W group_list=[YOUR PBS ACCOUNT GROUP ID]
#PBS -m e
#PBS -M user@email.com
```

### Compile GCHP
The first compile will take a long time, probably around 40 minutes give or take 20 minutes. On the head node, inside your run directory, run:
```
$ make build_all
```

This should end with 
```
###################################
###    GCHP executable exists!  ###
###################################
```
followed by some info about the run.

### Submit your job

Then submit:
```
$ qsub gchp.pleiades.run
```

The output netcdf files will be under OutputDir. I have configured the default output to include the v/v concentrations of the 3 CO2 "species" in the sample file I got from Seb and the meterological variables I thought would be necessary to interpret them. 

## Further information
I intend to make some of these steps less onerous in future, and this is still just a development branch of GCHP so things will change. I will do my best to keep this README up-to-date.

For further information about running GCHP, checkout the [wiki](http://wiki.seas.harvard.edu/geos-chem/index.php/Getting_Started_with_GCHP)

Also, feel free to email me at colin dot lee at dal dot ca
=======
