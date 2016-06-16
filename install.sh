#! /bin/bash

# Remove any prior installation
rm -rf ./bin

# Install the Kover python package and its dependencies
cd ./core
CORE=$(pwd)
(echo "Installing setuptools."; python bootstrap_setuptools.py; rm setuptools-*.zip)
cd kover/dataset/tools/contigs_tools
(echo "Building GATB tools"; sh install.sh)
cd ${CORE}
(echo "Installing Kover and its dependencies."; python setup.py install)
(echo "Cleaning up."; rm -r build; rm -r dist; rm -r kover.egg-info)
cd kover/dataset/tools/contigs_tools
rm -r build
cd ${CORE}
cd ..

# Create executable
mkdir bin
cd ./bin
(echo "Creating Kover executable."; ln -s ../interfaces/command_line.py kover; chmod u+x ./kover)
cd ..

echo Done.
