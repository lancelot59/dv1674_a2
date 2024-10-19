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


void correlation_coefficients(double*& result,std::vector<Vector>& datasets, double*& array, int dimension, int setstart, int setend)
{

    int count = (setstart * (2 * dimension - setstart - 1)) / 2; //calculate start position in result array

    for (int sample1 = setstart; sample1 < setend-1; sample1++) {
        for (auto sample2 = sample1 + 1 ; sample2 < dimension; sample2++) {
            result[count++] = pearson(datasets[sample1], datasets[sample2], array[sample1], array[sample2],array[sample1+dimension], array[sample2+dimension]);
        }
    }


}

double pearson(Vector vec1, Vector vec2, double x_mean, double y_mean, double x_mag, double y_mag) //this is multithreadable.
{ //performs relevant calculations with copied vectors
    vec1 - x_mean; 
    vec2 - y_mean;

    vec1 / x_mag ;
    vec2 / y_mag ;


    return std::max(std::min(vec1.dot(vec2), 1.0), -1.0);
}
};
