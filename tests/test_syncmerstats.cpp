#include "doctest.h"
#include "syncmerstats.hpp"
#include <stdexcept>

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
