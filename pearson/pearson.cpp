/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "analysis.hpp"
#include "dataset.hpp"
#include <iostream>
#include <cstdlib>

int main(int argc, char const* argv[])
{
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " [dataset] [outfile]" << std::endl;
        std::exit(1);
    }
    //we could if ive not missed anything entirely cut out std::vector for vector arrays. which should be way way quicker
    int dimension ;
    auto datasets { Dataset::read(argv[1],dimension) };
    std::vector<double> corrs = Analysis::correlation_coefficients(datasets, dimension);
    Dataset::write(corrs, argv[2]);

    return 0;
}
