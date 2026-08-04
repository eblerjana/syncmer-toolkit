#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <stdexcept>

extern "C" {
#include "syng.h"
}

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
}

int main(int argc, char* argv[]) {

	// parse command line
	string usage = "Usage: syncmer-toolkit <.1path file>";

	if (argc == 1) {
		// no arguments provided, just print usage info
		cerr << usage << endl;
		return 0;
	}

	if (argc < 2) {
		cerr << usage << "\n" << endl;
		cerr << "Error: Too few commandline arguments provided." << endl;
		return 1;
	}

	if (argc > 2) {
		cerr << usage << "\n" << endl;
		cerr << "Error: Too many commandline arguments provided." << endl;
		return 1;
	}

	string filename_path = argv[1];

	// make sure that files can be opened
	try {
		check_open_file(filename_path, "path");
	} catch (const runtime_error& e) {
		cerr << e.what();
		return 1;
	}

	cout << "Running program with the following files:" << endl;
	cout << "-----------------------------------------" << endl;
	cout << "1path file:\t" << filename_path << "\ņ" << endl;

	// compute syncmer stats
	int exit_code = compute_syncmer_stats(filename_path);

	return exit_code;
};
