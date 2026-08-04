#include <iostream>
#include <unordered_map>
#include "syncmerstats.hpp"

extern "C" {
#include "syng.h"
#include "seqhash.h"
}

using namespace std;

int compute_syncmer_stats (std::string& filename_path) {

	// read the 1path file
	OneSchema *schema = oneSchemaCreateFromText (syngSchemaText);
	OneFile* ipath = oneFileOpenRead(filename_path.data(), schema, "path", 1);

	if (!ipath) {
		cerr << "Error: could not open 1path file." << endl;
	}

	// maps each syncmer to its counts across all haplotypes
	vector<unordered_map<uint16_t, uint32_t>> counts;

	bool ok = oneReadLine(ipath);

	while(ok) {
		// read Z lines
		// iterate all nodes and increase counter for respective file ID
	}

	return 0;
}
