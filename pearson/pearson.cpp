/*
Author: David Holmqvist <daae19@student.bth.se>
*/
#include <thread>
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
    int dimension ; //=128
    auto datasets { Dataset::read(argv[1],dimension) };

    double* result = new double[dimension*dimension];
    int num_threads = 8;
    int setstart, setend;
    
    double* array = new double[dimension];//creating an array for all the mean values so they only need ot be calculated once
    for(auto sample = 0; sample < dimension; sample++ )
        array[sample] = datasets[sample].mean();

    std::vector<std::thread> threads;

    for(int thread_num = 0; thread_num < num_threads; thread_num++)
    {
        setstart = dimension*dimension/num_threads * thread_num; //= 2048 * threadnum
        setend = setstart + dimension*dimension/num_threads; //= 2048 *threadnum+1
        threads.emplace_back(Analysis::correlation_coefficients,std::ref(result),std::ref(datasets), std::ref(array), dimension,setstart,setend);
    }
    for(auto& thread: threads)
        thread.join();
    Dataset::write(result, argv[2], dimension);

    delete [] array;
    delete [] result;

    return 0;
}
