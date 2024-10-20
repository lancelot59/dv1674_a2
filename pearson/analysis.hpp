/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "vector.hpp"
#include <vector>

#if !defined(ANALYSIS_HPP)
#define ANALYSIS_HPP

namespace Analysis {
void correlation_coefficients(double*& result,std::vector<Vector>& datasets, int dimension, int setstart, int secondstart, int setend, int secondend);
};

#endif
