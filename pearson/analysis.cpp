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

std::vector<double> correlation_coefficients(std::vector<Vector> datasets)//this should be pretty easy to multithread
{
    std::vector<double> result {};

    for (auto sample1 { 0 }; sample1 < datasets.size() - 1; sample1++) {
        for (auto sample2 { sample1 + 1 }; sample2 < datasets.size(); sample2++) {
            auto corr { pearson(datasets[sample1], datasets[sample2]) };
            result.push_back(corr);
        }
    }

    return result;
}

double pearson(Vector vec1, Vector vec2) //this is probably multithreadable but im unsure if we get very much value from it we might lose out from the function call overheads.
{
    auto x_mean { vec1.mean() };
    auto y_mean { vec2.mean() };

    vec1 - x_mean;
    vec2 - y_mean;

    auto x_mag { vec1.magnitude() };
    auto y_mag { vec2.magnitude() };

    vec1 / x_mag ;
    vec2 / y_mag ;


    return std::max(std::min(vec1.dot(vec2), 1.0), -1.0);
}
};
