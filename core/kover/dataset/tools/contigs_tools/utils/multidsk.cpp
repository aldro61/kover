#include <DSK.hpp>
#include <DSK.cpp>
#include <fstream>

using namespace std;

/********************************************************************************/

int main (int argc, char* argv[])
{
	string line;
	string path(argv[2]);
	size_t found;
	
	// Searching for absolute path
	found = path.find_last_of("/");
	path = path.substr(0, found);
	
	// Opening file containing list of files to process
	ifstream genome_list (argv[2]);
	try
	{
		if (genome_list)
		{
			while (getline(genome_list, line))
			{
				// Path to file
				line = path + "/" + line;
				
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
