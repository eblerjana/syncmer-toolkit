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

void visit(int32_t& node, vector<int32_t>& counts, unordered_set<long long int>& seen) {
	long long int sync_id = abs(node);
	if (seen.find(sync_id) != seen.end()) {
		counts[sync_id] = -1;
	} else {
		increase_syncmer_count(counts, sync_id);
		seen.insert(sync_id);
	}
}

int compute_syncmer_stats_from_paths (string& pathfile_path, string& khashfile_path, string& outfile_path) {

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

int compute_syncmer_stats_from_gbwt (string& gbwtfile_path, string& khashfile_path, string& outfile_path) {
	
	// read the gbwt file
	OneSchema *schema = oneSchemaCreateFromText (syngSchemaText);
	OneFile* ofGBWT = oneFileOpenRead(gbwtfile_path.data(), schema, "gbwt", 1);
	
	if (!ofGBWT) {
		cerr << "Error: could not open gbwt file." << endl;
		oneSchemaDestroy(schema);
		return 1;
	}

	SyngBWT* sgb = syngBWTread(ofGBWT);
	oneFileClose(ofGBWT);

	if (!sgb) {
		cerr << "Error: coult not read GBWT from file " << gbwtfile_path << endl;
		oneSchemaDestroy(schema);
		return 1;
	}

	// reopen the file here again to be on the safe side (because we don't know what effect GBWT construction had on of GBWT)
	OneFile* ofPZ = oneFileOpenRead(gbwtfile_path.data(), schema, "gbwt", 1);

	if (!ofPZ) {
		cerr << "Error: could not re-open gbwt file." << endl;
		oneSchemaDestroy(schema);
		syngBWTdestroy(sgb);
		return 1;
	}

	// parse the P and Z lines to keep track of each paths source ID (which would otherwise be lost)
	struct PathStart {int32_t startNode; long long int j0; };
	vector<PathStart> starts;
	vector<int64_t> sourceFiles;
	int64_t source_id = 0;
	bool line_read = oneReadLine(ofPZ);

	while(line_read && ofPZ->lineType != 'V') {
		switch (ofPZ->lineType) {
			case 'P' : {
				long long int new_source_id = oneInt(ofPZ, 1);
		  		if (new_source_id != source_id) {
					if (new_source_id != source_id + 1) {
						cerr << "Error: paths in input file are not ordered." << endl;
						oneSchemaDestroy(schema);
						oneFileClose(ofPZ);
						syngBWTdestroy(sgb);
						return 1;
					}
					source_id = new_source_id;
				}
				break;
			} case 'Z': {
				long long int start_node = oneInt(ofPZ, 0);
				long long int j0 = oneInt(ofPZ, 2);
				PathStart path;
				path.startNode = start_node;
				path.j0 = j0;
				starts.push_back(path);
				sourceFiles.push_back(source_id);
				break;
			}
			default: break;
		}
		oneReadLine(ofPZ);
	}
	oneFileClose(ofPZ);
	oneSchemaDestroy(schema);

	// load the syncmers from khash file
	SyncmerSet *sms = syncmerSetRead(khashfile_path.data());
	int synLen = sms->params.w + sms->params.k;
	// max syncmer ID
	long long int nSyncmers = kmerHashMax (sms->kh);

	// keep a count for each syncmer. Entries of -1 indicate non-unique syncmers
	vector<int32_t> counts(nSyncmers + 1, 0);

	// nodes seen in a file so far
	unordered_set<long long int> seen;

	// traverse the paths through the GBWT
	int64_t current_file = -1;

	for (size_t i = 0; i < starts.size(); ++i) {
		if (sourceFiles[i] != current_file) {
			// entering next file
			seen.clear();
			current_file = sourceFiles[i];
		}

		// visit the start node of the path
		visit(starts[i].startNode, counts, seen);

		// traverse the rest of the path, starting from the start node
		SyngBWTpath *sbp = syngBWTpathStartOld(sgb, starts[i].startNode, starts[i].j0);

		int32_t nextNode; uint32_t offset;
		while (syngBWTpathNext(sbp, &nextNode, &offset)) {
			visit(nextNode, counts, seen);
		}
		syngBWTpathDestroy(sbp);
	}	
	syngBWTdestroy(sgb);

	// write out unique syncmers and their total counts	
	char *buf = new char[synLen + 1]();

	ofstream outfile;
	outfile.open(outfile_path + "_syncmers.tsv");
	if (!outfile.good()) {
		cerr << "Error: output file " << outfile_path << "_syncmers.tsv cannot be created." << endl;
		delete[] buf;
		syncmerSetDestroy(sms);
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
