# Install the Kover python package
cd ./src
python bootstrap_setuptools.py
python setup.py install
rm -r build
rm -r dist
rm -r kover.egg-info
cd ..

# Create executable
mkdir bin
cd ./bin
(echo "Creating Kover executable."; ln -s ../src/interfaces/command_line.py kover; chmod u+x ./kover)
cd ..

echo Done.
