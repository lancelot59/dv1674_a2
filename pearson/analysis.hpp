/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "vector.hpp"
#include <vector>

#if !defined(ANALYSIS_HPP)
#define ANALYSIS_HPP

namespace Analysis {
void correlation_coefficients(double*& result,std::vector<Vector>& datasets, double*& array, int dimension, int setstart, int setend);
};

#endif
