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
    int dimension;
    int setstart;
    int secondstart;
    int setend;
    int secondend;
    ThreadData()
        : result(nullptr), datasets(*(new std::vector<Vector>)), dimension(0), setstart(0), secondstart(0), setend(0), secondend(0) {}

    ThreadData(double* result, std::vector<Vector>& datasets,  int dimension,int setstart, int secondstart ,int setend, int secondend)
    :result(result),datasets(datasets),dimension(dimension),setstart(setstart),secondstart(secondstart),setend(setend), secondend(secondend)
    {}
    ThreadData(const ThreadData&) = delete;
    ThreadData& operator=(const ThreadData& input)
    {
        result = input.result;
        datasets = input.datasets;
        dimension = input.dimension;
        setstart = input.setstart;
        setend = input.setend;
        return *this;
    }
};

void* calculate_coefficients(void* arg)
{
    ThreadData* data = static_cast<ThreadData*>(arg);
    Analysis::correlation_coefficients(data->result, data->datasets, data->dimension, data->setstart, data->secondstart, data->setend, data->secondend);
    return nullptr;
}

int main(int argc, char const* argv[])
{
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " [dataset] [outfile] [numberofthreads]" << std::endl;
        std::exit(1);
    }
    //we could if ive not missed anything entirely cut out std::vector for vector arrays. which should be way way quicker
    int dimension ;
    std::vector<Vector> datasets { Dataset::read(argv[1],dimension) }; //read from file

    //until here 128 copy allocations


    int arraysize = dimension*(dimension-1)/2; //calculate array size
    double* result = new double[arraysize]; //create array


    int num_threads = std::stoi(argv[3]); //set number of threads
    int setstart, secondstart, setend, secondend;
    
    for(auto i = 0; i < dimension; i++ ) //perform precalculations
    {
        datasets[i] - datasets[i].mean();
        datasets[i] / datasets[i].magnitude();
    }

    pthread_t threads[num_threads];
    int blockSize = dimension/num_threads/2; //calculate block size (needs improved implementation)
    std::vector<ThreadData*> thread_data;

    for(int thread_num = 0; thread_num < num_threads; thread_num++)
    {
        setstart = blockSize * thread_num; //start of block for thread 
        setend = setstart + blockSize + 1;  //end of block for thread
        secondstart = dimension - blockSize - blockSize*thread_num;
        secondend = dimension - blockSize*thread_num; 
        if(thread_num != 0)
            secondend +=1;

        ThreadData* thisthread = new ThreadData{result, datasets,  dimension, setstart, secondstart, setend, secondend};
        thread_data.push_back(thisthread);
        if(pthread_create(&threads[thread_num], nullptr, calculate_coefficients,(void*)thisthread) != 0)
        {
            std::cerr << "error creating thread " << thread_num << std::endl;
            delete [] result;
            return 1;
        }
    }
    for(int i = 0; i < num_threads; i++) //terminate threads
        pthread_join(threads[i],nullptr);
    
    Dataset::write(result, argv[2], arraysize);


    for(auto& i:thread_data)
        delete i;
    delete [] result;

    return 0;
}
