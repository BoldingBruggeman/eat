# EAT
Ensemble and Assimilation Tool (EAT). 



The documentation is work in progress and is not complete. The instructions on how to download and compile should contain enough information to complete.



## Software installation

The EAT software is available from a GitHub repository. Please follow these instructions to download the core EAT code. In addition to this code the PDAF code must also be installed.

Here *~/source/repos/EAT* is used a installation folder as an example. Any folder can be used.

The commands can be copy pasted - with proper respect paid to the folders used.

```bash
mkdir ~/soure/repos/EAT && cd ~/soure/repos/EAT
git clone --recurse-submodules https://github.com/BoldingBruggeman/eat
```

Next step is to prepare the PDAF code-base to integrate with EAT.

This involves three steps.

Step 1: register and download the PDAF code

The get access to the PDAF source code [registration](http://pdaf.awi.de/register/index.php) is required. Please follow the instructions on the registration page and proceed to step 2 when the code is downloaded and installed.

Step 2: link the PDAF-code to the EAT code.

```bash
cd eat/extern
ln -s ~/PDAF-D_V1.16 pdaf
cd ..
```

The PDAF source code folder given in the above command must be changed to fit the actual folder.

Step 3: Copy PDAF CMake configuration files from EAT to PDAF.

PDAF does not come with CMake based configuration files. The following will prepare the PDAF source code to integrate and build together with EAT.

```bash
cp PDAF_CMake/CMakeLists.txt extern/pdaf
cp PDAF_CMake/src/CMakeLists.txt extern/pdaf/src
cp PDAF_CMake/src/pdaf_configure.h.in extern/pdaf/src
```

The above PDAF related steps are necessary for now - but might be changed/relaxed at a later stage.

## Compilation

CMake advocates 'out of source' compilation. The build-folder can be any folder. Here the simplest solution is used.

```bash
cd ../ && mkdir build && cd build
cmake ../eat
make install
```

An easy way to support multiple compilers on a single platform the build directory can include the compiler name like - *build/gfortran* or *build/ifort*. 

If the above completes without errors software configuration, compilation and installation is done successfully. An executable _eat_filter_pdaf_ has been installed into a user configurable folder (the default is platform dependent). 

## Running a test case

EAT comes with 2 test cases:

- [The SEAMLESS development case](tests/seamless). This case is used to develop and test the communication between  observation, filter and model components
- [The 2d case](tests/2d). This case resembles the standard 2D case from the [PDAF tutorials](http://pdaf.awi.de/files/pdaf_tutorial_the_model.pdf).

Another test case is provided with the [GOTM](https://www.gotm.net) EAT model component.

The cases are available in *.../tests* folder. Each of the tests has a README for specific information.

https://github.com/18F/open-source-guide/blob/18f-pages/pages/making-readmes-readable.md
