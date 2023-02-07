# EAT: Ensemble and Assimilation Tool

(*) The documentation is work in progress and is not complete. The instructions on how to download and compile should contain enough information to complete.

## Installation and use

Instructions for using EAT with GOTM-FABM are given on [the EAT wiki](https://github.com/BoldingBruggeman/eat/wiki/)

## Design

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

EAT comes with two fully functional configurations. 

1. The first - [2d](tests/2d) is a re-implementation of the example used in the [PDAF tutorials](http://pdaf.awi.de/files/pdaf_tutorial_the_model.pdf).
2. The second is uses GOTM as the model-component - with EAT interface code [here](models/gotm) and model configuration in [nns_annual](tests/nns_annual).
