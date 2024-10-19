/*
Author: David Holmqvist <daae19@student.bth.se>
*/
#include <thread>
#include "analysis.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <list>
#include <vector>

namespace Analysis {


void correlation_coefficients(double*& result,std::vector<Vector>& datasets, int dimension, int setstart, int setend)
{

    int count = (setstart * (2 * dimension - setstart - 1)) / 2; //calculate start position in result array

    for (int sample1 = setstart; sample1 < setend-1; sample1++) {
        for (auto sample2 = sample1 + 1 ; sample2 < dimension; sample2++) {
            result[count++] = std::max(std::min(datasets[sample1].dot(datasets[sample2]),1.0),-1.0);
        }
    }
}
};
