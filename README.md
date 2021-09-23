# EAT
Ensemble and Assimilation Tool (EAT). 



(*) The documentation is work in progress and is not complete. The instructions on how to download and compile should contain enough information to complete.

The basic idea of EAT is to ease the ensemble and data-assimilaton workflow by providing a framework separating the communication between components (observations, filter and model) from the actual filter-calculation, observation handling and model integration. The exchange between the different components is done exclusively via MPI.

The three components in EAT are:

1. the observation handler
2. the assimilation filter
3. the dynamical model

The responsibility of the observation handler is: 1) split the model integration in time slices according to when observations are available, 2) handle and pre-process observations to make them ready for use by the assimilation filter.

The assimilation component will - based on pre-processed observations and a list of forecasted state-vectors calculated analysed state-vectors.

The model component will do model integration in the time slices controlled by the observation handler and send state-vectors to the assimilation component.

The aim of EAT is that minimal - if any - changes are needed for different use cases and that e.g. the filter component can be used unchanged for a number of different models. EAT provides a thin wrapper layer around the full featured software components.

EAT does not enforce any specific programmig language and in the GOTM observation handler described below is implemented in Python.

As an example - the present implementation of the filter component uses [PDAF](http://pdaf.awi.de/trac/wiki) as the underlying library for the actual filter operations. The actual EAT wrapper only contains ~120 lines of Fortran code (excluding the required PDAF interface routines). The executable *eat_filter_pdaf* is 100% independent of both observation handler and integration model.

In principle other data-assimilation frameworks can be used as long as they implement the information protocol defined between the different components. For software written in Fortran this is easily achieved using the *eat_config.F90* module provided by EAT.

The EAT requirements for the the model component are that the model uses the - *initialize, integrate, finalize* workflow. Furthermore, the model must support different run-time configurations via configuration files and finally packing and un-packing of model fields and variables to/from a 1-dimensional state-vector must be implemented.

EAT is executed via a call to *mpiexec* with appropriate options as will be demonstrated later.

EAT contains minimal configuration but *eat_filter_pdaf* and Fortran programs using *eat_config.F90* can use *namelist-files* for configuration of e.g. output verbosity.

## Software installation

The EAT software is available from a GitHub repository. Please follow these instructions to download the core EAT code. In addition to this code the PDAF code must also be installed.

The command snippets given below can be copy pasted verbatim to a terminal. Windows users must adjust the workflow according the the Windows way.

For convenience 3 shell variables are introduced to hold different folders used during software cloning, software compilation and software insttalation:

`reposdir=~/source/repos/EAT`

Where the EAT source code is cloned from the GitHub repository.

`builddir=~/source/repos/EAT/build`

The folder where the actual compilation is done.

`bindir=~/local/eat/bin`

The folder where the compiled executables are installed. This folder can be added to the PATH environment variable.

The following commands will create the *$reposdir* folder and clone the EAT code - including any necessary Git submodules.

```bash
mkdir -p $reposdir && cd $reposdir
git clone --recurse-submodules https://github.com/BoldingBruggeman/eat
```

The next step is to prepare the [PDAF source code](http://pdaf.awi.de/trac/wiki) to be integrated with EAT.

Installation of PDAF involves three steps.

Step 1: register and download the PDAF code

The get access to the PDAF source code [registration](http://pdaf.awi.de/register/index.php) is required. Please follow the instructions given on the registration page and proceed to step 2 when the code is downloaded and installed.

Step 2: make the PDAF-code available to EAT

There are two methods to specify where the PDAF code is.

The first method is to create a link to point to the downloaded and un-packed PDAF source code:

```bash
cd $reposdir/eat/extern
ln -s ~/PDAF-D_V1.16 pdaf
cd ..
```

The PDAF source code folder given in the above link command must be adjusted to fit the actual folder where the PDAF code is downloaded.

Note the name of the linked folder must be *pdaf*.

The next method is to provide an additional argument to the *cmake* command as shown below.

Step 3: Copy PDAF CMake configuration files from EAT to PDAF.

PDAF does not come with CMake based configuration files. The following will prepare the PDAF source code to integrate and build with EAT.

```bash
cp PDAF_CMake/CMakeLists.txt extern/pdaf/
cp PDAF_CMake/src/CMakeLists.txt extern/pdaf/src/
cp PDAF_CMake/src/pdaf_configure.h.in extern/pdaf/src/
cp PDAF_CMAKE/src/pdaf_configure.F90 extern/pdaf/src/
```

Again - the folder where the files are copied must be adjusted to fit where the PDAF source code is installed. In the example above a link has been created to *$reposdir/eat/extern/pdaf*.

The above PDAF related steps are necessary for now - but might be changed/relaxed at a later stage.

## Compilation

EAT uses [CMake](https://www.cmake.org) to configure and compile the software. CMake advocates 'out of source' compilation. The build-folder can be any folder. Here the simplest solution is used.

```bash
cmake -B ../build
cmake --build ../build --target install
```

To specify the folder with the PDAF source code on the command line - use:

```bash
cmake -DEAT_PDAF_BASE=<folder_with_PDAF_code> -B ../build
cmake --build ../build --target install
```

The first command will configure the software and write build configuration files in the folder *../build*.

The default configuration of EAT includes building the GOTM EAT model component. For further information abut the GOTM component see [here](./models/gotm/).

The second command will build (compile) and install the software in the default installation folder. The default installation folder is  *~/local/eat/* - but it can be changed using the CMake configuration option *CMAKE_INSTALL_PREFIX*.

An easy way to support multiple compilers on a single platform the build directory can include the compiler name like - *build/gfortran* or *build/ifort*. 

If the above completes without errors software configuration, compilation and installation is done successfully. 

To use the installed executables the preferred way is to add e.g. *~/local/eat/bin* to the PATH environment variable. As an alternative (also used in this documentation) a variable *$bindir* is used.

## Included EAT ready to run cases

EAT comes with two fully functional configurations. 

1. The first - [2d](tests/2d) is a re-implementation of the example used in the [PDAF tutorials](http://pdaf.awi.de/files/pdaf_tutorial_the_model.pdf).
2. The second is uses GOTM as the model-component - with EAT interface code [here](models/gotm) and model configuration in [nns_annual](tests/nns_annual).
