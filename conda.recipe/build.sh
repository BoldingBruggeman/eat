#!/bin/bash

mkdir build
cd build

declare -a CMAKE_PLATFORM_FLAGS
if [[ ${HOST} =~ .*darwin.* ]]; then
  CMAKE_PLATFORM_FLAGS+=(-DCMAKE_OSX_SYSROOT="${CONDA_BUILD_SYSROOT}")
else
  CMAKE_PLATFORM_FLAGS+=(-DCMAKE_TOOLCHAIN_FILE="${RECIPE_DIR}/cross-linux.cmake")
fi

cmake $RECIPE_DIR/.. -DPython3_EXECUTABLE="$PYTHON" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$PREFIX ${CMAKE_PLATFORM_FLAGS[@]}
cmake --build . --config Release --parallel $CPU_COUNT --target install
cd ..
rm -r build
