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
 * Gael Letarte




EndOfMessage


# Remove any prior installation
rm -r ./bin 2>/dev/null

# Get the install context
INSTALL=$(pwd -P)
cd ./core && CORE=$(pwd -P)

echo -e "Installing GATB\n-----------------"
cd kover/dataset/tools/kmer_tools
sh install.sh
cd ${CORE}
echo -e "\n\n"

echo -e "Installing setuptools\n-----------------"
pip install setuptools
echo -e "\n\n"

echo -e "Installing Kover\n-----------------"
pip install .
echo -e "\n\n"

# Create executable
cd ..
mkdir bin
cd ./bin
(echo -e "Creating Kover executable"; ln -s ../interfaces/command_line.py kover; chmod u+x ./kover)
echo -e Done.
echo -e "\n\nPlease add $INSTALL/bin/ to your PATH environment variable."
