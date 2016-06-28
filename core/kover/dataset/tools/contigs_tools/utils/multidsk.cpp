/*
*	Kover: Learn interpretable computational phenotyping models from k-merized genomic data
*	Copyright (C) 2015  Alexandre Drouin & Gaël Letarte St-Pierre
*
*	This program is free software: you can redistribute it and/or modify
*	it under the terms of the GNU General Public License as published by
*	the Free Software Foundation, either version 3 of the License, or
*	(at your option) any later version.
*
*	This program is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*
*	You should have received a copy of the GNU General Public License
*	along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <DSK.hpp>
#include <DSK.cpp>
#include <fstream>

using namespace std;

/********************************************************************************/

int main (int argc, char* argv[])
{
	string line;
	
	// Opening file containing list of files to process
	ifstream genome_list (argv[2]);
	try
	{
		if (genome_list)
		{
			while (getline(genome_list, line))
			{
				
				// Modifying path argv for DSK
				argv[2] = strdup(line.c_str());
				
				// Applying DSK tool to the file
				DSK().run(argc, argv);
			}
		}
	}
	catch (Exception& e)
	{
		genome_list.close();
		cerr << "EXCEPTION: " << e.getMessage() << endl;
		return EXIT_FAILURE;
	}
	// Closing file
	genome_list.close();
    return EXIT_SUCCESS;
}
