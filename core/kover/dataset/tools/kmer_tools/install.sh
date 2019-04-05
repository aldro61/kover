# Remove any prior build
rm -rf build 2>/dev/null

# Preparing gatb-core folder
cd thirdparty
rm -rf gatb-core 2>/dev/null
mkdir gatb-core
cd gatb-core

#Download gatb-core
echo '* Downloading gatb-core'
case "$(python --version 2>&1)" in
    *" 2."*)
        python -c 'import urllib; urllib.urlretrieve ("https://github.com/GATB/gatb-core/archive/v1.2.2.tar.gz", "gatbcore.tar.gz")'
        ;;
	*" 3."*)
		python -c 'import urllib.request; urllib.request.urlretrieve ("https://github.com/GATB/gatb-core/archive/v1.2.2.tar.gz", "gatbcore.tar.gz")'
	;;
    *)
        echo "Error: Python 2 or 3 is not installed"
        ;;
esac
tar -zxf gatbcore.tar.gz --strip-components=1
rm gatbcore.tar.gz
cd ../..
# Create build directory
mkdir build
cd build

# Prepare the makefile
echo '* Compiling gatb-core'
cmake .. -Wno-dev

# Run the makefile to build dsk, multidsk and dsk2kover tools:
make -j8

# Clean up gatb-core source files
echo '* Cleaning up gatb-core source files'
cd ..
cd thirdparty
rm -rf gatb-core
cd ..

echo "* Copying binaries to their target location"
mv build/bin/* .
