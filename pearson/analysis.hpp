/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "vector.hpp"
#include <vector>

#if !defined(ANALYSIS_HPP)
#define ANALYSIS_HPP

namespace Analysis {
void correlation_coefficients(double*& result,std::vector<Vector>& datasets, double*& array, int dimension, int setstart, int setend);
double pearson(Vector vec1, Vector vec2, double x_mean, double y_mean, double x_mag, double y_mag);
};

#endif
