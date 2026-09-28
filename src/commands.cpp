#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include "syncmerstats.hpp"
#include "commands.hpp"
#include "syng_wrapper.hpp"

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


int create_histogram(int argc, char* argv[]) {

	string usage = "Usage:\tsyncmer-toolkit histogram <.1path file> <.1khash file> <outname>\n\tsyncmer-toolkit histogram <.1gbwt file> <.1khash file> <outname>\n";
	if (argc < 5) {
		cerr << usage << "\n" << endl;
		cerr << "Error: Too few commandline arguments provided." << endl;
		return 1;
	}

	if (argc > 5) {
		cerr << usage << "\n" << endl;
		cerr << "Error: Too many commandline arguments provided." << endl;
		return 1;
	}

	string sourcefile_path = argv[2];
	string khashfile_path = argv[3];
	string outfile_path = argv[4];

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

	cerr << "Running program with the following files:" << endl;
	cerr << "-----------------------------------------" << endl;
	cerr << ftype << " file:\t" << sourcefile_path << endl;
	cerr << "1khash file:\t" << khashfile_path << endl << endl;

	int exit_code = 0;

	// compute syncmer stats
	if (ftype == "path") {
		exit_code = compute_syncmer_stats_from_paths(sourcefile_path, khashfile_path, outfile_path);
	} else {
		exit_code = compute_syncmer_stats_from_gbwt(sourcefile_path, khashfile_path, outfile_path);
	}

	return exit_code;
};

int compute_distances(int argc, char* argv[]) {

	string usage = "Usage:\tsyncmer-toolkit distances <.1path file> <.1khash file> <syncmerset> <outname>\n";

	if (argc < 6) {
		cerr << usage << endl;
		cerr << "Error: Too few arguments provided." << endl;
		return 1;
	}

	if (argc > 6) {
		cerr << usage << endl;
		cerr << "Error: Too many arguments provided." << endl;
		return 1;
	}

	string sourcefile_path = argv[2];
	string khashfile_path = argv[3];
	string syncmerset_path = argv[4];
	string outfile_path = argv[5];

	// make sure that 1path and 1khash files can be opened
	try {
		check_open_file(sourcefile_path, "path");
		check_open_file(khashfile_path, "khash");
	} catch (const runtime_error& e) {
		cerr << e.what();
		return 1;
	}

	cerr << "Running program with the following files:" << endl;
	cerr << "-----------------------------------------" << endl;
	cerr << "path file:\t" << sourcefile_path << endl;
	cerr << "khash file:\t" << khashfile_path << endl;
	cerr << "syncmerset file:\t" << syncmerset_path << endl << endl;

	// compute distances
	return compute_syncmer_distances_from_paths(sourcefile_path, khashfile_path, syncmerset_path, outfile_path);
}
