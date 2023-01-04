# Linux

`conda-build <RECIPEDIR> -c conda-forge`

Adding `--override-channels --no-test` may be necessary to work around an MPI
conflict reported during the testing phase. This happened 2022-11-29 with Python
3.9 targets; Python 3.8 was fine. Update 2023-01-04: this seems to have been fixed.

# Windows

`conda-build <RECIPEDIR> -c conda-forge`

Before this command, the VS 2017 environment needs to be loaded.
To do this, use "x64 Native Tools Command Prompt for VS 2017" in the start menu.