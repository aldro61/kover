#! /bin/bash


cat << EndOfMessage

██╗  ██╗ ██╗██╗  ██╗   ██╗███████╗██████╗ 
██║ ██╔╝██╔╝╚██╗ ██║   ██║██╔════╝██╔══██╗
█████╔╝██╔╝  ╚██╗██║   ██║█████╗  ██████╔╝
██╔═██╗╚██╗  ██╔╝╚██╗ ██╔╝██╔══╝  ██╔══██╗
██║  ██╗╚██╗██╔╝  ╚████╔╝ ███████╗██║  ██║
╚═╝  ╚═╝ ╚═╝╚═╝    ╚═══╝  ╚══════╝╚═╝  ╚═╝
                                          
By: Alexandre Drouin

Contributors:
 * Gaël Letarte Saint-Pierre                                      




EndOfMessage


# Remove any prior installation
rm -r ./bin 2>/dev/null

cd ./core
CORE=$(pwd -P)

echo -e "Installing GATB\n-----------------"
cd kover/dataset/tools/kmer_tools
sh install.sh
cd ${CORE}
echo -e "\n\n"

echo -e "Installing setuptools\n-----------------"
python bootstrap_setuptools.py; rm setuptools-*.zip
echo -e "\n\n"

echo -e "Installing Kover\n-----------------"
python setup.py install
echo -e "\n\n"

echo -e "Cleaning up..."
cd kover/dataset/tools/kmer_tools
rm -r build
cd ${CORE}

# Create executable
mkdir bin
cd ./bin
(echo -e "Creating Kover executable"; ln -s ../interfaces/command_line.py kover; chmod u+x ./kover)
cd ..

echo -e Done.
