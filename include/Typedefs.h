#ifndef Typedefs_header
#define Typedefs_header

#include <string>
#include <vector>
#include <array>
#include <chrono>

typedef std::chrono::high_resolution_clock::time_point timePoint;
typedef std::vector<std::string> strVec_t;
typedef std::vector<std::vector<std::string>> strVecVec_t;
typedef std::vector<std::array<std::string, 2>> strArrVec_t;
typedef std::array<std::string, 2> strArr_t;
typedef std::array<std::vector<std::string>, 2> strVecArr_t;
typedef std::array<std::vector<std::vector<std::string>>, 2> strVecVecArr_t;

#endif
