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

    int arraysize = dimension*(dimension-1)/2;
    double* result = new double[arraysize];
    int num_threads = dimension;
    int setstart, setend;
    
    double* array = new double[dimension];//creating an array for all the mean values so they only need ot be calculated once
    for(auto sample = 0; sample < dimension; sample++ )
        array[sample] = datasets[sample].mean();

    std::vector<std::thread> threads;

    int blockSize = dimension/num_threads;

    for(int thread_num = 0; thread_num < num_threads; thread_num++)
    {
        setstart = blockSize * thread_num; //start of block for thread 
        setend = setstart + blockSize;  //end of block for thread    
        if(thread_num != num_threads-1)
            setend += 1;

        threads.emplace_back(Analysis::correlation_coefficients,std::ref(result),std::ref(datasets), std::ref(array), dimension,setstart,setend);
    }
    for(auto& thread: threads) //terminate threads
        thread.join();
/*
    int unused = 0;
    for(int i = 0; i < dimension*dimension; i++)
    {
        if(result[i] != 0 && unused != 0)
            result[unused++] = result[i];
        else if( result[i] == 0 && unused == 0)
            unused = i;
    }
*/    
    
    Dataset::write(result, argv[2], arraysize);



    delete [] array;
    delete [] result;

    return 0;
}
