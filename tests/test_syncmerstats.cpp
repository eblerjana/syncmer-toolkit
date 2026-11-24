#include "doctest.h"
#include "syncmerstats.hpp"
#include <stdexcept>
#include "test_helpers.hpp"
#include <map>
#include <sstream>
#include <iostream>

using namespace std;

TEST_CASE("Test syncmerstats::count: counting.") {
	vector<int32_t> counts = {0,0,0};
	increase_syncmer_count(counts, 1);
	CHECK(counts[1] == 1);
	counts[0] = -1;
	increase_syncmer_count(counts, 0);
	CHECK(counts[0] == -1);
}

TEST_CASE("Test syncmerstats::visit: counting.") {
	vector<int32_t> counts = {0,2,3,6};
	unordered_set<long long int> seen;
	int32_t node = 2;
	visit(node, counts, seen);
	CHECK(counts[2] == 4);
	CHECK(seen.find(2) != seen.end());
}

TEST_CASE("Test syncmerstats::visit: mark syncmer as non-unique.") {
	vector<int32_t> counts(5,0);
	unordered_set<long long int> seen;
	int32_t first = 4;
	int32_t second = -4;
	visit(first, counts, seen);
	CHECK(counts[4] == 1);
	CHECK(seen.find(4) != seen.end());
	visit(second, counts, seen);
	CHECK(counts[4] == -1);
	CHECK(seen.find(-4) == seen.end());
}

TEST_CASE("Test syncmerstats::get_max: read sample file.") {
	string filename = string(TEST_DATA_DIR) + "small.1khash";
	CHECK(get_max(filename) == 11);
}

TEST_CASE("Test syncmerstats::get_max missing file.") {
	string non_existent = "non-existent.1khash";
	CHECK_THROWS_AS(get_max(non_existent), runtime_error);
}


TEST_CASE("Test syncmerstats::compute_syncmer_anchors_from_paths ") {
	string path = data_path("small.1path");
	string subset = data_path("small_subset.tsv");

	// easy example with two haplotypes
	string path_example1 = data_path("anchors1.1path");
	string subset_example1 = data_path("anchors1_subset.tsv");

	// multiple haplotypes
	string path_example2 = data_path("anchors2.1path");
	string subset_example2 = data_path("anchors2_subset.tsv");

	// reference path split to multiple contigs
	string path_example3 = data_path("anchors3.1path");
	string subset_example3 = data_path("anchors3_subset.tsv");

	// some haplotypes are reverse complemented
	string path_example4 = data_path("anchors4.1path");
	string subset_example4 = data_path("anchors4_subset.tsv");

	// some anchors are reversed
	string path_example5 = data_path("anchors5.1path");
	string subset_example5 = data_path("anchors5_subset.tsv");

	// haplotypes split in contigs, some with no candidates
	string path_example6 = data_path("anchors6.1path");
	string subset_example6 = data_path("anchors6_subset.tsv");

	// reference paths are not listed first in file, named paths
	string path_example7 = data_path("anchors7.1path");
	string subset_example7 = data_path("anchors7_subset.tsv");

	// no haplotypes, just reference path in graph
	string path_example8 = data_path("anchors8.1path");
	string subset_example8 = data_path("anchors8_subset.tsv");


	/* test if issues with missing input files are handled properly
	 */

	SUBCASE("non-existent reference sequence name") {
		string out = output_path("non-existent.tsv");
		string reference = "non-existant";
		REQUIRE(compute_syncmer_anchors_from_paths(path, subset, reference, out) == 1);
	}

	SUBCASE("output file cannot be created") {
		string out = output_path("non/existent/dir/output.tsv");
		string reference = "seq.fa";
		REQUIRE(compute_syncmer_anchors_from_paths(path, subset, reference, out) == 1);
	}

	SUBCASE("input path file cannot be opened") {
		string out = output_path("output.tsv");
		string non_existent_path = "non-existent.1path";
		string reference = "seq.fa";
		REQUIRE(compute_syncmer_anchors_from_paths(non_existent_path, subset, reference, out) == 1);
	}

	SUBCASE("subset file cannot be opened") {
		string out = output_path("output.tsv");
		string reference = "seq.fa";
		string non_existent_subset = "non-existent.tsv";
		REQUIRE(compute_syncmer_anchors_from_paths(path, non_existent_subset, reference, out) == 1);	
	}

	/*
	 * test basic functionality of the anchors algorithm, with corner cases covered.
	 * */

	SUBCASE("test anchors example 1") {
		string out = output_path("anchors1.tsv");
		string reference = "ref.fa";
		REQUIRE(compute_syncmer_anchors_from_paths(path_example1, subset_example1, reference, out) == 0);
		string expected_output = 
			"syncmer_ID\treference_contig\treference_position\n"
			"1\t1\t2\n"
			"7\t1\t62\n"
			"8\t1\t72\n"
			"12\t1\t112\n";
		REQUIRE(read_file(out) == expected_output);
	}

	SUBCASE("test anchors example 2") {
		string out = output_path("anchors2.tsv");
		string reference = "ref.fa";
		REQUIRE(compute_syncmer_anchors_from_paths(path_example2, subset_example2, reference, out) == 0);
		string expected_output = 
			"syncmer_ID\treference_contig\treference_position\n"
			"1\t1\t2\n"
			"2\t1\t12\n"
			"4\t1\t32\n"
			"5\t1\t42\n"
			"7\t1\t62\n"
			"12\t1\t112\n";
		REQUIRE(read_file(out) == expected_output);
	}

	SUBCASE("test anchors example 3") {
		string out = output_path("anchors3.tsv");
		string reference = "ref.fa";
		REQUIRE(compute_syncmer_anchors_from_paths(path_example3, subset_example3, reference, out) == 0);
		string expected_output = 
			"syncmer_ID\treference_contig\treference_position\n"
			"6\t1\t22\n"
			"7\t1\t32\n";
		REQUIRE(read_file(out) == expected_output);
	}

	SUBCASE("test anchors example 4") {
		string out = output_path("anchors4.tsv");
		string reference = "ref.fa";
		REQUIRE(compute_syncmer_anchors_from_paths(path_example4, subset_example4, reference, out) == 0);
		string expected_output = 
			"syncmer_ID\treference_contig\treference_position\n"
			"1\t1\t2\n"
			"2\t1\t12\n"
			"5\t1\t42\n"
			"6\t1\t52\n"
			"7\t1\t62\n";
		REQUIRE(read_file(out) == expected_output);
	}

	SUBCASE("test anchors example 5") {
		string out = output_path("anchors5.tsv");
		string reference = "ref.fa";
		REQUIRE(compute_syncmer_anchors_from_paths(path_example5, subset_example5, reference, out) == 0);
		string expected_output = 
			"syncmer_ID\treference_contig\treference_position\n"
			"1\t1\t2\n"
			"-2\t1\t12\n"
			"4\t1\t32\n"
			"7\t1\t62\n"
			"8\t1\t72\n"
			"9\t2\t2\n"
			"10\t2\t12\n"
			"11\t2\t22\n";
		REQUIRE(read_file(out) == expected_output);
	}

	SUBCASE("test anchors example 6") {
		string out = output_path("anchors6.tsv");
		string reference = "ref.fa";
		REQUIRE(compute_syncmer_anchors_from_paths(path_example6, subset_example6, reference, out) == 0);
		string expected_output = 
			"syncmer_ID\treference_contig\treference_position\n"
			"1\t1\t2\n"
			"2\t1\t12\n"
			"3\t1\t22\n"
			"4\t2\t2\n"
			"5\t2\t12\n"
			"6\t2\t22\n";
		REQUIRE(read_file(out) == expected_output);
	}

	SUBCASE("test anchors example 7") {
		string out = output_path("anchors7.tsv");
		string reference = "ref.fa";
		REQUIRE(compute_syncmer_anchors_from_paths(path_example7, subset_example7, reference, out) == 0);
		string expected_output = 
			"syncmer_ID\treference_contig\treference_position\n"
			"1\tchrB\t0\n"
			"2\tchrB\t5\n"
			"3\tchrB\t17\n"
			"4\t1\t3\n"
			"5\t1\t10\n"
			"6\t1\t24\n"
			"7\t1\t31\n";
		REQUIRE(read_file(out) == expected_output);
	}

	SUBCASE("test anchors example 8") {
		string out = output_path("anchors8.tsv");
		string reference = "ref.fa";
		REQUIRE(compute_syncmer_anchors_from_paths(path_example8, subset_example8, reference, out) == 0);
		string expected_output = 
			"syncmer_ID\treference_contig\treference_position\n"
			"1\t1\t2\n"
			"3\t1\t22\n"
			"4\t1\t42\n";
		REQUIRE(read_file(out) == expected_output);
	}
}


// =============================================================================
// The following test cases were generated by Claude (see: tests/data/README.md)
// =============================================================================

static const string SMALL_SYNCMERS =
	"syncmer_ID\ttotal_count\n"
	"1\t0\n" "2\t0\n" "3\t1\n" "4\t1\n" "5\t1\n" "6\t1\n"
	"7\t1\n" "8\t1\n" "9\t1\n" "10\t1\n" "11\t1\n";
static const string SMALL_HISTOGRAM = "0\t2\n" "1\t9\n" "2\t0\n";

static const string REPEATS_SYNCMERS =
	"syncmer_ID\ttotal_count\n"
	"1\t0\n" "2\t0\n" "3\t2\n" "4\t1\n" "5\t1\n" "6\t1\n"
	"7\t0\n" "9\t0\n" "10\t1\n" "11\t1\n";
static const string REPEATS_HISTOGRAM = "0\t4\n" "1\t5\n" "2\t1\n";

static const string SMALL_DISTANCES =
	"source_id\tdistance\tcount\n"
	"1\t13\t1\n" "1\t11\t1\n" "1\t4\t1\n" "1\t2\t1\n"
	"2\t12\t1\n" "2\t11\t1\n" "2\t2\t4\n" "2\t1\t1\n";

static const string SMALL_DISTANCES_EMPTY_SUBSET =
	"source_id\tdistance\tcount\n"
	"1\t17\t1\n" "1\t13\t1\n"
	"2\t17\t1\n" "2\t15\t1\n";

static const string START_AT_ZERO_DISTANCES =
	"source_id\tdistance\tcount\n"
	"1\t13\t1\n" "1\t11\t1\n" "1\t6\t1\n"
	"2\t12\t1\n" "2\t11\t1\n" "2\t2\t4\n" "2\t1\t1\n";

static map<long long, long long> total_length_per_source(const string& filename) {
	map<long long, long long> totals;
	istringstream in(read_file(filename));
	string header;
	getline(in, header);
	long long source, distance, count;
	while (in >> source >> distance >> count) totals[source] += distance * count;
	return totals;
}

TEST_CASE("Test syncmerstats::get_max: khash file without t lines.") {
	string filename = data_path("small_no_t_lines.1khash");
	CHECK_THROWS_WITH_AS(get_max(filename), "Error: khash file does not contain any t lines.", runtime_error);
}

TEST_CASE("Test syncmerstats::compute_syncmer_stats_from_paths") {
	string path = data_path("small.1path");
	string khash = data_path("small.1khash");

	SUBCASE("sample files") {
		string out = output_path("stats_paths");
		REQUIRE(compute_syncmer_stats_from_paths(path, khash, out) == 0);
		CHECK(read_file(out + "_syncmers.tsv") == SMALL_SYNCMERS);
		CHECK(read_file(out + "_histogram.tsv") == SMALL_HISTOGRAM);
	}

	SUBCASE("shared and non-unique syncmers") {
		string repeats = data_path("small_repeats.1path");
		string out = output_path("stats_paths_repeats");
		REQUIRE(compute_syncmer_stats_from_paths(repeats, khash, out) == 0);
		CHECK(read_file(out + "_syncmers.tsv") == REPEATS_SYNCMERS);
		CHECK(read_file(out + "_histogram.tsv") == REPEATS_HISTOGRAM);
	}

	SUBCASE("missing path file") {
		string missing = data_path("non-existent.1path");
		string out = output_path("stats_paths_missing_path");
		CHECK(compute_syncmer_stats_from_paths(missing, khash, out) == 1);
	}

	SUBCASE("missing khash file") {
		string missing = data_path("non-existent.1khash");
		string out = output_path("stats_paths_missing_khash");
		CHECK(compute_syncmer_stats_from_paths(path, missing, out) == 1);
	}

	SUBCASE("paths not ordered by source file") {
		string unordered = data_path("small_unordered.1path");
		string out = output_path("stats_paths_unordered");
		CHECK(compute_syncmer_stats_from_paths(unordered, khash, out) == 1);
	}

	SUBCASE("output file cannot be created") {
		string out = output_path("non-existent-dir/stats");
		CHECK(compute_syncmer_stats_from_paths(path, khash, out) == 1);
	}
}

TEST_CASE("Test syncmerstats::compute_syncmer_stats_from_gbwt") {
	string gbwt = data_path("small.1gbwt");
	string khash = data_path("small.1khash");

	SUBCASE("sample files give the same result as the path version") {
		string out = output_path("stats_gbwt");
		REQUIRE(compute_syncmer_stats_from_gbwt(gbwt, khash, out) == 0);
		CHECK(read_file(out + "_syncmers.tsv") == SMALL_SYNCMERS);
		CHECK(read_file(out + "_histogram.tsv") == SMALL_HISTOGRAM);
	}

	SUBCASE("missing gbwt file") {
		string missing = data_path("non-existent.1gbwt");
		string out = output_path("stats_gbwt_missing");
		CHECK(compute_syncmer_stats_from_gbwt(missing, khash, out) == 1);
	}

	SUBCASE("missing khash file") {
		string missing = data_path("non-existent.1khash");
		string out = output_path("stats_gbwt_missing_khash");
		CHECK(compute_syncmer_stats_from_gbwt(gbwt, missing, out) == 1);
	}

	SUBCASE("output file cannot be created") {
		string out = output_path("non-existent-dir/stats");
		CHECK(compute_syncmer_stats_from_gbwt(gbwt, khash, out) == 1);
	}
}

TEST_CASE("Test syncmerstats::compute_syncmer_distances_from_paths") {
	string path = data_path("small.1path");
	string khash = data_path("small.1khash");
	string subset = data_path("small_subset.tsv");

	SUBCASE("sample files") {
		string out = output_path("distances.tsv");
		REQUIRE(compute_syncmer_distances_from_paths(path, khash, subset, out) == 0);
		CHECK(read_file(out) == SMALL_DISTANCES);
	}

	SUBCASE("pieces of each source file add up to its total path length") {
		string out = output_path("distances_totals.tsv");
		REQUIRE(compute_syncmer_distances_from_paths(path, khash, subset, out) == 0);
		map<long long, long long> totals = total_length_per_source(out);
		CHECK(totals.size() == 2);
		CHECK(totals[1] == 17 + 13);
		CHECK(totals[2] == 17 + 15);
	}

	SUBCASE("empty syncmer set: every path is a single piece") {
		string empty = data_path("small_subset_empty.tsv");
		string out = output_path("distances_empty_set.tsv");
		REQUIRE(compute_syncmer_distances_from_paths(path, khash, empty, out) == 0);
		CHECK(read_file(out) == SMALL_DISTANCES_EMPTY_SUBSET);
	}

	SUBCASE("no piece before a syncmer at position 0") {
		string zero = data_path("small_start_at_zero.1path");
		string out = output_path("distances_start_at_zero.tsv");
		REQUIRE(compute_syncmer_distances_from_paths(zero, khash, subset, out) == 0);
		CHECK(read_file(out) == START_AT_ZERO_DISTANCES);
	}

	SUBCASE("paths not ordered by source file") {
		string unordered = data_path("small_unordered.1path");
		string out = output_path("distances_unordered.tsv");
		CHECK(compute_syncmer_distances_from_paths(unordered, khash, subset, out) == 1);
	}

	SUBCASE("o line without preceding P line") {
		string broken = data_path("small_o_without_p.1path");
		string out = output_path("distances_o_without_p.tsv");
		CHECK(compute_syncmer_distances_from_paths(broken, khash, subset, out) == 1);
	}

	SUBCASE("missing path file") {
		string missing = data_path("non-existent.1path");
		string out = output_path("distances_missing_path.tsv");
		CHECK(compute_syncmer_distances_from_paths(missing, khash, subset, out) == 1);
	}

	SUBCASE("missing khash file") {
		string missing = data_path("non-existent.1khash");
		string out = output_path("distances_missing_khash.tsv");
		CHECK(compute_syncmer_distances_from_paths(path, missing, subset, out) == 1);
	}

	SUBCASE("missing syncmer set file") {
		string missing = data_path("non-existent.tsv");
		string out = output_path("distances_missing_set.tsv");
		CHECK(compute_syncmer_distances_from_paths(path, khash, missing, out) == 1);
	}

	SUBCASE("syncmer ID larger than number of syncmers") {
		string invalid = data_path("small_subset_invalid.tsv");
		string out = output_path("distances_invalid_set.tsv");
		CHECK(compute_syncmer_distances_from_paths(path, khash, invalid, out) == 1);
	}

	SUBCASE("o line without preceding z line") {
		string broken = data_path("small_o_without_z.1path");
		string out = output_path("distances_o_without_z.tsv");
		CHECK(compute_syncmer_distances_from_paths(broken, khash, subset, out) == 1);
	}

	SUBCASE("output file cannot be created") {
		string out = output_path("non-existent-dir/distances.tsv");
		CHECK(compute_syncmer_distances_from_paths(path, khash, subset, out) == 1);
	}
}

// small.1path has no I lines, so its paths are addressed by their 1-based "sequence number in
// file" (the third field of the P line). small_with_names.1path is the same data but with an I
// line (sequence name) added to each path - see tests/data/README.md. Both fixtures share the
// same syncmer content, so expected outputs below apply to whichever fixture matches the region.
//
//  seq.fa  (source 1), path 1/"chr1" (len 17): z = 3, -4, -5 ; o = 2, 4, 6
//  seq.fa  (source 1), path 2/"chr2" (len 13): z = -6        ; o = 2
//  seq2.fa (source 2), path 1/"chr1" (len 17): z = -7, -8, -9; o = 2, 4, 6
//  seq2.fa (source 2), path 2/"chr2" (len 15): z = -10, 11   ; o = 2, 3
//
// small_subset.tsv contains syncmers 3, 5, 7, 8, 9, 10, 11 (all with count 1); syncmers 4, 6
// are not in the subset even though they lie on these paths.
TEST_CASE("Test syncmerstats::compute_syncmer_query_from_paths") {
	string path = data_path("small.1path");
	string named = data_path("small_with_names.1path");
	string khash = data_path("small.1khash");
	string subset = data_path("small_subset.tsv");
	string seqfa = "seq.fa";
	string seq2fa = "seq2.fa";

	SUBCASE("index-based addressing on a file with no I lines; forward and reverse orientation") {
		// covers: happy path via index fallback, forward-traversed syncmer 3, reverse-traversed syncmer 5
		string out = output_path("query_index_based.tsv");
		string region = "1:0-17";
		REQUIRE(compute_syncmer_query_from_paths(path, khash, subset, seqfa, region, out) == 0);
		CHECK(read_file(out) ==
			"syncmer_ID\tposition\tcount\tsequence\n"
			"3\t2\t1\tggctcgcta\n"
			"-5\t6\t1\tcgctacgat\n");
	}

	SUBCASE("name-based addressing on a file with I lines; all syncmers reverse-traversed") {
		string out = output_path("query_named.tsv");
		string region = "chr1:0-17";
		REQUIRE(compute_syncmer_query_from_paths(named, khash, subset, seq2fa, region, out) == 0);
		CHECK(read_file(out) ==
			"syncmer_ID\tposition\tcount\tsequence\n"
			"-7\t2\t1\tgcctcgctt\n"
			"-8\t4\t1\tctcgcttcg\n"
			"-9\t6\t1\tcgcttcgat\n");
	}

	SUBCASE("the same sequence name on a different source file resolves to a different path") {
		string out = output_path("query_named_other_source.tsv");
		string region = "chr1:0-17";
		REQUIRE(compute_syncmer_query_from_paths(named, khash, subset, seqfa, region, out) == 0);
		CHECK(read_file(out) ==
			"syncmer_ID\tposition\tcount\tsequence\n"
			"3\t2\t1\tggctcgcta\n"
			"-5\t6\t1\tcgctacgat\n");
	}

	SUBCASE("region before any syncmer on the path: header-only output") {
		string out = output_path("query_before.tsv");
		string region = "1:0-1";
		REQUIRE(compute_syncmer_query_from_paths(path, khash, subset, seqfa, region, out) == 0);
		CHECK(read_file(out) == "syncmer_ID\tposition\tcount\tsequence\n");
	}

	SUBCASE("region after all syncmers on the path: header-only output") {
		string out = output_path("query_after.tsv");
		string region = "1:7-17";
		REQUIRE(compute_syncmer_query_from_paths(path, khash, subset, seqfa, region, out) == 0);
		CHECK(read_file(out) == "syncmer_ID\tposition\tcount\tsequence\n");
	}

	SUBCASE("syncmer inside the region but absent from the subset file is excluded") {
		// syncmer 4 sits at position 4 on path 1, inside [3,5), but is not in small_subset.tsv
		string out = output_path("query_not_in_subset.tsv");
		string region = "1:3-5";
		REQUIRE(compute_syncmer_query_from_paths(path, khash, subset, seqfa, region, out) == 0);
		CHECK(read_file(out) == "syncmer_ID\tposition\tcount\tsequence\n");
	}

	SUBCASE("unknown fasta filename") {
		string out = output_path("query_unknown_fasta.tsv");
		string region = "1:0-17";
		string unknown_fasta = "doesnotexist.fa";
		CHECK(compute_syncmer_query_from_paths(path, khash, subset, unknown_fasta, region, out) == 1);
	}

	SUBCASE("unknown sequence index on a file with no I lines") {
		string out = output_path("query_unknown_index.tsv");
		string region = "99:0-17";
		CHECK(compute_syncmer_query_from_paths(path, khash, subset, seqfa, region, out) == 1);
	}

	SUBCASE("unknown sequence name on a file with I lines") {
		string out = output_path("query_unknown_name.tsv");
		string region = "chrX:0-17";
		CHECK(compute_syncmer_query_from_paths(named, khash, subset, seqfa, region, out) == 1);
	}

	SUBCASE("non-numeric seqid against a file with no I lines") {
		string out = output_path("query_non_numeric_seqid.tsv");
		string region = "chr1:0-17";
		CHECK(compute_syncmer_query_from_paths(path, khash, subset, seqfa, region, out) == 1);
	}

	SUBCASE("malformed region: missing colon") {
		string out = output_path("query_malformed_colon.tsv");
		string region = "0-17";
		CHECK(compute_syncmer_query_from_paths(path, khash, subset, seqfa, region, out) == 1);
	}

	SUBCASE("malformed region: non-integer start") {
		string out = output_path("query_malformed_start.tsv");
		string region = "1:abc-17";
		CHECK(compute_syncmer_query_from_paths(path, khash, subset, seqfa, region, out) == 1);
	}

	SUBCASE("malformed region: start >= end") {
		string out = output_path("query_malformed_order.tsv");
		string region = "1:17-5";
		CHECK(compute_syncmer_query_from_paths(path, khash, subset, seqfa, region, out) == 1);
	}

	SUBCASE("missing path file") {
		string missing = data_path("non-existent.1path");
		string out = output_path("query_missing_path.tsv");
		string region = "1:0-17";
		CHECK(compute_syncmer_query_from_paths(missing, khash, subset, seqfa, region, out) == 1);
	}

	SUBCASE("missing khash file") {
		string missing = data_path("non-existent.1khash");
		string out = output_path("query_missing_khash.tsv");
		string region = "1:0-17";
		CHECK(compute_syncmer_query_from_paths(path, missing, subset, seqfa, region, out) == 1);
	}

	SUBCASE("missing syncmer set file") {
		string missing = data_path("non-existent.tsv");
		string out = output_path("query_missing_set.tsv");
		string region = "1:0-17";
		CHECK(compute_syncmer_query_from_paths(path, khash, missing, seqfa, region, out) == 1);
	}

	SUBCASE("output file cannot be created") {
		string out = output_path("non-existent-dir/query.tsv");
		string region = "1:0-17";
		CHECK(compute_syncmer_query_from_paths(path, khash, subset, seqfa, region, out) == 1);
	}
}
