# Linux

`conda-build <RECIPEDIR> -c conda-forge`

Adding `--override-channels --no-test` may be necessary to work around an MPI
conflict reported during the testing phase. This happened 2022-11-29 with Python
3.9 targets; Python 3.8 was fine. Update 2023-01-04: this seems to have been fixed.

# Windows

`conda-build <RECIPEDIR> -c conda-forge`

If VS2017 is installed, this command should work on a normal command prompt.
(no need to load MSVC or Intel-specific environments)