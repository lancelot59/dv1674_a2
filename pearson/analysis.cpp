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
#include <mutex>

namespace Analysis {

std::mutex result_mutex;

void correlation_coefficients(double*& result,std::vector<Vector>& datasets, double*& array, int dimension, int setstart, int setend)
{
    int count = setstart; //This is wrong

    for (int sample1 = setstart; sample1 < setend-1; sample1++) {
        for (auto sample2 = sample1 + 1 ; sample2 < dimension; sample2++) {
            result[count++] = pearson(datasets[sample1], datasets[sample2], array[sample1], array[sample2]);
        }
    }


}

double pearson(Vector vec1, Vector vec2, double x_mean, double y_mean) //this is probably multithreadable but im unsure if we get very much value from it we might lose out from the function call overheads.
{
    vec1 - x_mean;
    vec2 - y_mean;

    x_mean = vec1.magnitude() ;
    y_mean =  vec2.magnitude() ;

    vec1 / x_mean ;
    vec2 / y_mean ;


    return std::max(std::min(vec1.dot(vec2), 1.0), -1.0);
}
};
