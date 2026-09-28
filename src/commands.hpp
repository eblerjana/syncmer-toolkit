#ifndef COMMANDS_HPP
#define COMMANDS_HPP

#include <string>

bool ends_with(std::string const & value, std::string const & ending);
void check_open_file(std::string filename, std::string filetype);
int create_histogram(int argc, char* argv[]);
int compute_distances(int argc, char* argv[]);

#endif // COMMANDS_HPP
