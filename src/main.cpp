#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include "syng_wrapper.hpp"

#include "syncmerstats.hpp"

using namespace std;

bool ends_with(string const & value, string const & ending) {
	/**
	 * Check if string (value) ends with string (ending)
	 * **/
	if (ending.size() > value.size()) return false;
	return equal(ending.rbegin(), ending.rend(), value.rbegin());
}


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
	string usage = "Usage:\tsyncmer-toolkit <.1path file> <.1khash file> <outname>\n\tsyncmer-toolkit <.1gbwt file> <.1khash file> <outname>\n";

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

	string sourcefile_path = argv[1];
	string khashfile_path = argv[2];
	string outfile_path = argv[3];

	string ftype = "";
	if (ends_with(sourcefile_path, ".1path")) {
		ftype = "path";
	} else if (ends_with(sourcefile_path, ".1gbwt")) {
		ftype = "gbwt";
	} else {
		cerr << "Error: first file must be either .1path or .1gbwt file." << endl;
		return 1;
	}

	// make sure that 1path file can be opened
	try {
		check_open_file(sourcefile_path, ftype);
		check_open_file(khashfile_path, "khash");
	} catch (const runtime_error& e) {
		cerr << e.what();
		return 1;
	}

	cout << "Running program with the following files:" << endl;
	cout << "-----------------------------------------" << endl;
	cout << ftype << " file:\t" << sourcefile_path << endl;
	cout << "1khash file:\t" << khashfile_path << endl << endl;

	int exit_code = 0;

	// compute syncmer stats
	if (ftype == "path") {
		exit_code = compute_syncmer_stats_from_paths(sourcefile_path, khashfile_path, outfile_path);
	} else {
		exit_code = compute_syncmer_stats_from_gbwt(sourcefile_path, khashfile_path, outfile_path);
	}

	return exit_code;
};
