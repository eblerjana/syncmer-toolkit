#ifndef SYNCMER_STATS_HPP
#define SYNCMER_STATS_HPP

#include <string>

int compute_syncmer_stats_from_paths (std::string& pathfile_path, std::string& khashfile_path, std::string& outfile_path);
int compute_syncmer_stats_from_gbwt (std::string& gbwtfile_path, std::string& khashfile_path, std::string& outfile_path);

#endif // SYNCMER_STATS_HPP
