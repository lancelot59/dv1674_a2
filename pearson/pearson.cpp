/*
Author: David Holmqvist <daae19@student.bth.se>
*/
#include <pthread.h>
#include "analysis.hpp"
#include "dataset.hpp"
#include <iostream>
#include <cstdlib>

struct ThreadData {
    double* result;
    std::vector<Vector>& datasets;
    double* array;
    int dimension;
    int setstart;
    int setend;
    ThreadData()
        : result(nullptr), datasets(*(new std::vector<Vector>)), array(nullptr), dimension(0), setstart(0), setend(0) {}

    ThreadData(double* result, std::vector<Vector>& datasets, double* array, int dimension,int setstart,int setend)
    :result(result),datasets(datasets),array(array),dimension(dimension),setstart(setstart),setend(setend)
    {}
    ThreadData(const ThreadData&) = delete;
    ThreadData& operator=(const ThreadData& input)
    {
        result = input.result;
        datasets = input.datasets;
        array = input.array;
        dimension = input.dimension;
        setstart = input.setstart;
        setend = input.setend;
        return *this;
    }
};

void* calculate_coefficients(void* arg)
{
    ThreadData* data = static_cast<ThreadData*>(arg);
    Analysis::correlation_coefficients(data->result, data->datasets, data->array, data->dimension, data->setstart, data->setend);
    return nullptr;
}

int main(int argc, char const* argv[])
{
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " [dataset] [outfile]" << std::endl;
        std::exit(1);
    }
    //we could if ive not missed anything entirely cut out std::vector for vector arrays. which should be way way quicker
    int dimension ;
    std::vector<Vector> datasets { Dataset::read(argv[1],dimension) }; //read from file

    int arraysize = dimension*(dimension-1)/2; //calculate array size
    double* result = new double[arraysize]; //create array
    int num_threads = 1; //set number of threads
    int setstart, setend;
    
    double* array = new double[dimension*2];//creating an array for all the mean values so they only need ot be calculated once
    for(auto sample = 0; sample < dimension; sample++ )
        array[sample] = datasets[sample].mean();
    for(auto sample = 0; sample < dimension; sample++)
        array[sample+dimension] = datasets[sample].magnitude();

    pthread_t threads[num_threads];
    ThreadData* thread_data = new ThreadData[num_threads];
    int blockSize = dimension/num_threads; //calculate block size (needs improved implementation)
    for(int thread_num = 0; thread_num < num_threads; thread_num++)
    {
        setstart = blockSize * thread_num; //start of block for thread 
        setend = setstart + blockSize;  //end of block for thread    
        if(thread_num != num_threads-1) //required for proper alignment
            setend += 1;
        thread_data[thread_num] = {result, datasets, array, dimension, setstart, setend};
        if(pthread_create(&threads[thread_num], nullptr, calculate_coefficients,(void*)&thread_data[thread_num]) != 0)
        {
            std::cerr << "error creating thread " << thread_num << std::endl;
            delete [] array;
            delete [] result;
            return 1;
        }
    }
    for(int i = 0; i < num_threads; i++) //terminate threads
        pthread_join(threads[i],nullptr);
    
    Dataset::write(result, argv[2], arraysize);



    delete [] array;
    delete [] result;

    return 0;
}
