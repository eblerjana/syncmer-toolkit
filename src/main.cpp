#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include "commands.hpp"
#include "syng_wrapper.hpp"

using namespace std;

int main(int argc, char* argv[]) {

	// parse command line
	string usage = "Usage:\tsyncmer-toolkit histogram <.1path file> <.1khash file> <outname>\n\tsyncmer-toolkit histogram <.1gbwt file> <.1khash file> <outname>\n\tsyncmer-toolkit distances <.1path file> <.1khash file> <syncmerset> <outname>\n";

	if (argc < 3) {
		// no arguments provided, just print usage info
		cerr << usage << endl;
		return 0;
	}

	string program = argv[1];

	if (program == "histogram") {
		return create_histogram(argc, argv);
	} else if (program == "distances") {
		return compute_distances(argc, argv);
	} else {
		cerr << usage << "\n" << endl;
		cerr << "Error: program name must be specified as histogram or distances."  << endl;
		return 1;
	}
};
