### A simple 2D example

As an example and as proof of concept the '2d' case from the [PDAF documentation](http://pdaf.awi.de/files/pdaf_tutorial_the_model.pdf) has be implemented in EAT. The implementation is very similar to the original PDAF case but with the following differences:

- All files are generated - i.e. there are no data files provided by the set-up.
- The observation points are randomly distributed in the domain and the number of points varies between the _nobsmin_ and _nobsmax_  variables in the _eat_2d.nml_ namelist file.
- The standard deviation of the observations can be specified in the *eat_2d.nml* namelist file.

The model simulation is a progressing wave through a rectangular basin. The difference between the ensemble members is all in the initial conditions where a phase shift is introduced. 

![eat_2d](/home/kb/source/repos/EAT/eat/tests/2d/eat_2d_true.gif)

Observations are generated from the true solution by adding normal distributed random noise.

#### Compiling the 2d case:

All commands will be executed in the folder of the test.

```
cd $reposdir/eat/tests/2d
```

It is not strictly necessary to have built and installed  EAT - as it will be done as part of the command give below (if not already done). 

The compilation of the test case is done as follows:

```bash
cmake --build ../../../build --target test_eat_2d install
```

A successful compilation will have generated 3 executables - _eat_pdaf_filter_, _eat_obs_2d_ and _eat_model_2d_. The *eat_filter_pdaf* have been installed but *eat_obs_2d* and *eat_model_2d*  will be available in the folders where their respective source files are in the build folder.

#### Running the 2d case

The following commands assumes that OpenMPI has been used. As _mpiexec_ is not standardized using a different implementation of MPI will likely require adjustments to command examples given.

First - to get the true solution - run the following command:

```bash
mpiexec -np 1 $builddir/tests/2d/eat_obs_2d : -np 1 $builddir/tests/2d/eat_model_2d
```

Note that both the observation and model programs must be executed.

This will produce 19 _true\_??????????????.dat_ files and 17 _obs\_??????????????.dat_ files. Note - as the observation files are generated randomly this will **not** be the observations used by the assimilation done below.

To do the full assimilation the command below must be executed:

```bash
mpiexec --oversubscribe -np 1 $builddir/tests/2d/eat_obs_2d : -np 1 $bindir/eat_filter_pdaf : -np 9 $builddir/tests/2d/eat_model_2d
```

Note that above command is one long line and it must be copied to the terminal verbatim. The number of ensemble members **must** be 9 in accordance with the implementation from PDAF. 

Note also that _--oversubscribed_ is required by OpenMPI to start more processes than cores on the computer (here 11).

The command will use _mpiexec_ to start 3 programs - the _filter_-program, the  _observation_-program and the _model_-program. The first 2 each start 1 process and the last start 9 processes.

The command will write output from all processes to the screen and it can be difficult to trace progress. There are two methods available that can increase the clarity. 

The first method is to set the namelist variable _all\_verbose=.false._ in the _eat\_2d.nml_ namelist file (it defaults to true). Then only the first of the model processes will generate output. 

The second method is to use OpenMPI's version of _mpiexec_'s command line option - _--output-filename logs_. In this case the output are written in the folder _logs_ sorted according to process number.

```bash
mpiexec --oversubscribe --output-filename logs -np 1 $builddir/tests/2d/eat_obs_2d : -np 1 $bindir/eat_filter_pdaf : -np 9 $builddir/tests/2d/eat_model_2d
```

The command will generate a large number of files - for each time of observations  forecasted and analyzed fields for each ensemble member. The naming scheme of the files contains the type of file (*true, obs, for* and *ana* ) as well a a time stamp.

#### Plot/animate the 2d case results

A Python program has been included for creating plots and animations of the generated 2d-fields.

```bash
kb@orca:~/source/repos/EAT/eat/tests/2d$ python plot_fields.py -h
usage: plot_fields.py [-h] [--noshow] [--save] [--animate] [--interval INTERVAL] [--title TITLE] [infiles [infiles ...]]

Plot/animate 2D field(s).

positional arguments:
  infiles

optional arguments:
  -h, --help           show this help message and exit
  --noshow             do not show on screen
  --save               save png/mp4 file(s)
  --animate            animate input files
  --interval INTERVAL  interval (ms) between updates
  --title TITLE        animation title

```

To view the true solution on the screen:

```bash
./plot_fields.py true_??????????????.dat
```

To view the observation fields on the screen:

```bash
./plot_fields.py obs_??????????????.dat
```

To save .pngs of the fields:

```bash
./plot_fields.py obs_??????????????.dat --noview --save
```

To create an online animation of the fields:

```bash
./plot_fields.py obs_????????????.dat --animate
```

Note it is also possible to save to a file - either animated GIF or MP4 (note this might fail due to missing Python requirements):
```bash
./plot_fields.py obs_??????????????.dat --animate --save
```

![eat_2d_obs](/home/kb/source/repos/EAT/eat/tests/2d/eat_2d_obs.gif)
