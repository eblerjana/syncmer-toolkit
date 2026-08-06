#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>
#include <unordered_set>
#include "syncmerstats.hpp"

#include "syng_wrapper.hpp"

using namespace std;

void increase_syncmer_count(vector<int32_t>& counts, size_t position) {
	if (counts[position] < 0) {
		// non unique, don't add.
		return;
	} else {
		counts[position] += 1;
	}	
}

int compute_syncmer_stats (string& pathfile_path, string& khashfile_path, string& outfile_path) {

	// read the 1path file
	OneSchema *schema = oneSchemaCreateFromText (syngSchemaText);
	OneFile* ipath = oneFileOpenRead(pathfile_path.data(), schema, "path", 1);

	if (!ipath) {
		cerr << "Error: could not open 1path file." << endl;
		oneSchemaDestroy(schema);
		oneFileClose(ipath);
		return 1;
	}

	// read the khash file to look up syncmer counts and  sequences
	// use this as a template: https://github.com/richarddurbin/syng/blob/main/syngmap.c
	SyncmerSet *sms = syncmerSetRead(khashfile_path.data());
	int synLen = sms->params.w + sms->params.k;
	// max syncmer ID
	long long int nSyncmers = kmerHashMax (sms->kh);

	// keep a count for each syncmer. Entries of -1 indicate non-unique syncmers
	vector<int32_t> counts(nSyncmers + 1, 0);

	// use this as a template: https://github.com/richarddurbin/syng/blob/main/syngpath2gbwt.c
	bool line_read = oneReadLine(ipath);
	int64_t source_id = 0;
	// nodes seen in a file so far
	unordered_set<long long int> seen;

	while(line_read) {
		switch (ipath->lineType) {
			case 'P': {
				// get source file index from P lines
				int64_t new_source_id = oneInt(ipath, 1);
				if (new_source_id != source_id) {
					if (source_id % 50 == 0) cout << "Start reading paths from file " << new_source_id << " ..." << endl;
					if (new_source_id != source_id + 1) {
						cout << source_id << " " << new_source_id << endl;
						cerr << "Error: paths in input file are not ordered." << endl;
						oneSchemaDestroy(schema);
						oneFileClose(ipath);
						return 1;
					}
					// source file changed, therefore reset set of seen elements
					seen.clear(); 
					source_id = new_source_id;
				}
				break;
			}
			case 'z': {
				// from z lines, get list of syncmer nodes traversed by path
				long long int n_sync = oneLen(ipath);
				long long int* syncs = oneIntList(ipath);
				// iterate all syncmer nodes
				for (size_t i = 0; i < n_sync; ++i) {
					long long int sync_id = std::abs(syncs[i]);
					if (seen.find(sync_id) != seen.end()) {
						// set global count to -1, syncmer is not unique
						counts[sync_id] = -1;
					} else {
						// update global count
						increase_syncmer_count(counts, sync_id);
						// mark node as seen
						seen.insert(sync_id);
					}
				}
				break;
			}
			default:
				break;
		}
		line_read = oneReadLine(ipath);
	}
	oneFileClose(ipath);

	// write out unique syncmers and their total counts	
	char *buf = new char[synLen + 1]();

	ofstream outfile;
	outfile.open(outfile_path + "_syncmers.tsv");
	if (!outfile.good()) {
		cerr << "Error: output file " << outfile_path << "_syncmers.tsv cannot be created." << endl;
		delete[] buf;
		oneSchemaDestroy(schema);
		oneFileClose(ipath);
		return 1;
	}

	// write header and prepare lines for syncmers not covered by any file
	size_t total_unique = 0;
	outfile << "syncmer_ID\tsyncmer_seq_canonical\ttotal_count" << endl;

	// record how often each count (0,1,...,maxfile_ID) was seen
	vector<uint32_t> histogram(source_id + 1, 0);

	// write counts for each syncmer
	for (size_t sync_id = 1; sync_id < counts.size(); ++sync_id) {
		// do not output syncmers that are not unique
		if (counts[sync_id] < 0) continue;
		// look up syncmer sequence
		total_unique += 1;
		char* seq = kmerHashSeq(sms->kh, sync_id, buf);
		outfile << sync_id << "\t" << seq << "\t" << counts[sync_id] << endl;
		assert (counts[sync_id] < histogram.size());		
		histogram[counts[sync_id]] += 1;
	}
	outfile.close();
	delete[] buf;
	syncmerSetDestroy(sms);
	oneSchemaDestroy(schema);

	// write out the histogram file as well
	outfile.open(outfile_path + "_histogram.tsv");
 	if (!outfile.good()) {
		cerr << "Error: output file " << outfile_path << "_histogram.tsv cannot be created." << endl;
		return 1;
	}
	for (size_t i = 0; i < histogram.size(); ++i) {
		outfile << i << "\t" << histogram[i] << endl;
	}
	outfile.close();

	cout << "Wrote syncmer statistics to " << outfile_path << endl;
	cout << "Total syncmers:\t" << nSyncmers << endl;
	cout << "Total unique syncmers:\t" << total_unique << endl;
	return 0;
}
