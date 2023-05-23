#!/bin/bash

mkdir build
cd build

external_model_dir=$RECIPE_DIR/../extern

default_institutes="akvaplan;au;bb;csiro;examples;gotm;iow;jrc;msi;nersc;niva;pclake;pml;selma;su;uhh"
external_institutes="ersem;ogs;pisces"
fabm_institutes="-DFABM_INSTITUTES=$default_institutes;$external_institutes"

model_dirs="-DFABM_ERSEM_BASE=$external_model_dir/ersem -DFABM_OGS_BASE=$external_model_dir/ogs -DFABM_PISCES_BASE=$external_model_dir/pisces"

declare -a CMAKE_PLATFORM_FLAGS
if [[ ${HOST} =~ .*darwin.* ]]; then
  CMAKE_PLATFORM_FLAGS+=(-DCMAKE_OSX_SYSROOT="${CONDA_BUILD_SYSROOT}")
else
  CMAKE_PLATFORM_FLAGS+=(-DCMAKE_TOOLCHAIN_FILE="${RECIPE_DIR}/cross-linux.cmake")
fi

cmake $RECIPE_DIR/.. -DPython3_EXECUTABLE="$PYTHON" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$PREFIX ${CMAKE_PLATFORM_FLAGS[@]} $fabm_institutes $model_dirs
cmake --build . --config Release --parallel $CPU_COUNT --target install
cd ..
rm -r build
