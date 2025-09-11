This directory contains all c++ code and libraries

# Builing the lib 
1. create build directory and do cd build, then from there:
cmake -DCMAKE_BUILD_TYPE=Debug -DBUILD_TESTING=ON ..
cmake --build . -- -j$(nproc)
ctest --output-on-failure
