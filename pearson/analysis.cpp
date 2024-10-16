/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "analysis.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <list>
#include <vector>

namespace Analysis {

std::vector<double> correlation_coefficients(std::vector<Vector>& datasets, int dimension)//this should be pretty easy to multithread
{
    std::vector<double> result {};
    double x_mean;
    double y_mean;
    auto i = datasets.size();


    for (auto sample1 { 0 }; sample1 < i - 1; sample1++) {
        for (auto sample2 { sample1 + 1 }; sample2 < i; sample2++) {
            auto corr { pearson(datasets[sample1], datasets[sample2], x_mean, y_mean) };
            result.push_back(corr);
        }
    }

    return result;
}

double pearson(Vector vec1, Vector vec2, double x_mean, double y_mean) //this is probably multithreadable but im unsure if we get very much value from it we might lose out from the function call overheads.
{
    x_mean = vec1.mean() ;
    y_mean = vec2.mean() ;

    vec1 - x_mean;
    vec2 - y_mean;

    x_mean = vec1.magnitude() ;
    y_mean =  vec2.magnitude() ;

    vec1 / x_mean ;
    vec2 / y_mean ;


    return std::max(std::min(vec1.dot(vec2), 1.0), -1.0);
}
};
