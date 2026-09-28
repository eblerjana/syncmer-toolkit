#ifndef SYNCMER_STATS_HPP
#define SYNCMER_STATS_HPP

#include <vector>
#include <unordered_set>
#include <string>

void increase_syncmer_count(std::vector<int32_t>& counts, size_t position);
void visit(int32_t& node, std::vector<int32_t>& counts, std::unordered_set<long long int>& seen);
long long int get_max(std::string& khash_filename);

int compute_syncmer_stats_from_paths (std::string& pathfile_path, std::string& khashfile_path, std::string& outfile_path);
int compute_syncmer_stats_from_gbwt (std::string& gbwtfile_path, std::string& khashfile_path, std::string& outfile_path);
int compute_syncmer_distances_from_paths (std::string& pathfile_path, std::string& khashfile_path, std::string& syncmerset_path, std::string& outfile_path);



#endif // SYNCMER_STATS_HPP
