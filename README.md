# eat
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

Next step is to prepare the PDAF code base to integrate with EAT.

This involves two steps.

Step 1: link the PDAF-code to the EAT code.

```bash
cd eat/extern
ln -s ~/PDAF-D_V1.16 pdaf
cd ..
```

The PDAF source code folder given in the above command must be changed to fit the actual folder.

Step 2: Copy PDAF CMake configuration files from EAT to PDAF.

```bash
cp PDAF_CMake/CMakeLists.txt extern/pdaf
cp PDAF_CMake/src/CMakeLists.txt extern/pdaf/src
cp PDAF_CMake/src/pdaf_configure.h.in extern/pdaf/src
```

The above steps are necessary for now - but might changed/relaxed at a later stage.

## Compilation

CMake advocates 'out of source' compilation. The build-folder can be any folder. Here the simplest solution is used.

```bash
cd ../ && mkdir build && cd build
cmake ../eat
make
```

If the above completes software installation, configuration and compilation is done successfully.

## Running a test case

EAT comes with 3 test cases:

- The SEAMLESS development case. This case is used to develop and test the communication between  filter, observation and model components
- The 2d case. This case resembles the standard 2D case from the PDAF tutorials.
- A GOTM case. This case is used as an example of how to create a data-assimilation setup for GOTM.

The cases are available in *.../tests* folder. Each of the tests has a README for specific information.
