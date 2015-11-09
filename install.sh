#! /bin/bash

# Install the Kover python package and its dependencies
cd ./src
(echo "Installing setuptools."; python bootstrap_setuptools.py; rm setuptools-*.zip)
(echo "Installing Kover and its dependencies."; python setup.py install)
(echo "Cleaning up."; rm -r build; rm -r dist; rm -r kover.egg-info)
cd ..

# Create executable
mkdir bin
cd ./bin
(echo "Creating Kover executable."; ln -s ../src/interfaces/command_line.py kover; chmod u+x ./kover)
cd ..

echo Done.
