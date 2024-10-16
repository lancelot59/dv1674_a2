/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "vector.hpp"
#include <string>
#include <vector>

#if !defined(DATASET_HPP)
#define DATASET_HPP

namespace Dataset
{
    std::vector<Vector> read(std::string filename, int& dimension);
    void write(double*& data, std::string filename, int dimension);
};

#endif
