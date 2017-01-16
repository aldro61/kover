# Remove any prior build
rm -rf build

# Create build directory
mkdir  build
cd build

# Prepare the makefile
cmake ..

# Run the makefile to build dsk, multidsk and dsk2kover tools:
make -j8
