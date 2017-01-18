# Remove any prior build
rm -rf build

#Removing gatb-core previous version
cd thirdparty
rm -rf gatb-core
mkdir gatb-core
cd gatb-core

#Downloading gatb-core
VERSION=v1.2.2
URL=https://github.com/GATB/gatb-core/archive/$VERSION.tar.gz
wget $URL
tar -zxf $VERSION.tar.gz --strip-components=1
rm $VERSION.tar.gz
cd ../..

# Create build directory
mkdir  build
cd build

# Prepare the makefile
cmake ..

# Run the makefile to build dsk, multidsk and dsk2kover tools:
make -j8
