#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include "syng_wrapper.hpp"

#include "syncmerstats.hpp"

using namespace std;

void check_open_file(string filename, string filetype) {
	/**
	 * Check if ONE file can be opened.
	 **/

	OneSchema *schema = oneSchemaCreateFromText(syngSchemaText);
	OneFile *ofK = oneFileOpenRead(filename.data(), schema, filetype.data(), 1);
	if (!ofK) {
		stringstream ss;
		ss << "Error: the file " << filename << " cannot be opened." << endl;
		throw runtime_error(ss.str());
	}
	oneFileClose(ofK);
	oneSchemaDestroy(schema);
}

int main(int argc, char* argv[]) {

	// parse command line
	string usage = "Usage: syncmer-toolkit <.1path file> <.1khash file> <outfile>";

	if (argc == 1) {
		// no arguments provided, just print usage info
		cerr << usage << endl;
		return 0;
	}

	if (argc < 4) {
		cerr << usage << "\n" << endl;
		cerr << "Error: Too few commandline arguments provided." << endl;
		return 1;
	}

	if (argc > 4) {
		cerr << usage << "\n" << endl;
		cerr << "Error: Too many commandline arguments provided." << endl;
		return 1;
	}

	string pathfile_path = argv[1];
	string khashfile_path = argv[2];
	string outfile_path = argv[3];

	// make sure that 1path file can be opened
	try {
		check_open_file(pathfile_path, "path");
		check_open_file(khashfile_path, "khash");
	} catch (const runtime_error& e) {
		cerr << e.what();
		return 1;
	}

	cout << "Running program with the following files:" << endl;
	cout << "-----------------------------------------" << endl;
	cout << "1path file:\t" << pathfile_path << endl;
	cout << "1khash file:\t" << khashfile_path << endl << endl;

	// compute syncmer stats
	int exit_code = compute_syncmer_stats(pathfile_path, khashfile_path, outfile_path);

	return exit_code;
};
