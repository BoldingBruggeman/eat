**NOTE - the eat_pdaf_filter is not ready yet so the actual assimilation command given below will fail**

### A simple 2D example

This test case resembles the '2d' case from the [PDAF documentation](http://pdaf.awi.de/files/pdaf_tutorial_the_model.pdf) withe the following differences:

- All files are generated - i.e. there are no data-files provided by the set-up.
- The observation points are randomly distributed on the domain and the number of points varies between the _nobsmin_ and _nobsmax_ in the _eat_2d.nml_ namelist file.

#### Compiling the 2d case

It is assumed that a general configuration and compilation of _eat_ has been done. 

The _builddir_ below is the user provided folder used by _CMake_ to store the necessary compilation configuration. 

Then the compilation of the test case is done as follows:

```bash
cd $builddir
make test_eat_2d
```

A successful compilation will have compiled 3 executables - _eat_pdaf_filter_, _eat_2d_obs_ and _eat_2d_model_. For the test cases there is no _make install_ and the exectutables will be in the folders where their respective source files are.

#### Running the 2d case

The _sourcedir_ is the folder where the _eat_ source code is installed via the _git clone_ -command.

The commands will be executed in the folder of the test.

```
cd $installdir/tests/2d
```

The following commands assumes that OpenMPI has been used. As _mpiexec_ is not standardized using a different implementation of MPI will likely require adjustments to command examples given.

First - to get the true solution - run the following command:

```
cd $installdir/tests/2d
mpiexec -np 1 $builddir/tests/2d/eat_2d_obs : -np 1 $builddir/tests/2d/eat_2d_model
```

Note that both the observation and model programs must be executed.

This will produce 19 _true\_????.dat_ files and 17 _obs\_????.dat_ files. Note - as the observation files are generated randomly this will **not** be the observations used by the assimilation done below.

To do the full assimilation the command below must be executed:

```bash
mpiexec --oversubscribe -np 1 $builddir/src/eat_pdaf_filter : -np 1 $builddir/tests/2d/eat_2d_obs : -np 9 $builddir/tests/2d/eat_2d_model
```

Note that above command is one long line and it must be copied to the terminal verbatim.

Note also that _--oversubscribed_ is required by OpenMPI to start more processes than cores on the computer.

The command will use _mpiexec_ to start 3 programs - a _filter_-program, an _observation_-program and a _model_-program. The first 2 each start 1 process and the last start 9 processes (in accordance with the configuration of the PDAF-setup).

The command will write output from all processes to the screen and it can be difficult to trace progress. There are two methods available that can increase the clarity. 

The first method is to set the namelist variable _all\_verbose=.false._ in the _eat\_2d.nml_ namelist file (it defaults to true). Then only the first of the model processes will generate output. 

The second method is to use OpenMPI's version of _mpiexec_'s command line option - _--output-filename logs_. In this case the output are written in the folder _logs_ sorted according to process number.

```bash
mpiexec --oversubscribe --output-filename logs -np 1 $builddir/src/eat_pdaf_filter : -np 1 $builddir/tests/2d/eat_2d_obs : -np 9 $builddir/tests/2d/eat_2d_model
```

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
./plot_fields.py true_????.dat
```

To view the observation fields on the screen:

```bash
./plot_fields.py obs_????.dat
```

To save .pngs of the fields:

```bash
./plot_fields.py obs_????.dat --noview --save
```

To create an online animation of the fields:

```bash
./plot_fields.py obs_????.dat --animate
```

Note it is also possible to save to a .mp4 file (note this might fail due to missing Python requirements):
```bash
./plot_fields.py obs_????.dat --animate --save
```

