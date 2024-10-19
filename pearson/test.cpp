#include <iostream>
#include "vector.hpp"

int main()
{
    double* array = new double[8];
    for(int i = 0; i < 8; i++)
        array[i] = i+10;
    double* array2 = new double[8];
    for(int i = 0; i < 8; i++)
        array2[i] = i+10;
    Vector vec1 = {8,array}; 
    Vector vec2 = {8,array}; 
    int mean = vec1.mean();
    int mag1 = vec1.magnitude();

    vec1 - mean;
    vec2 - mean;
     
     int mag2 = vec2.magnitude();

     vec1/mag1;
     vec2/mag2;


     std::cout << "vec1: " << vec1.dot(vec1) << std::endl;
     std::cout << "vec2: " << vec2.dot(vec2) << std::endl;

     delete [] array;
     
}