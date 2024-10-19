/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "vector.hpp"
#include <iostream>
#include <cmath>
#include <vector>
#include <immintrin.h>

Vector::Vector()
    : size{0}, data{nullptr}
{
}

Vector::~Vector()
{
    if (data)
    {
        delete[] data;
    }

    size = 0;
}

Vector::Vector(unsigned size)
    : size{size}, data{new double[size]}
{
}

Vector::Vector(unsigned size, double *data)
    : size{size}, data{data}
{
}

Vector::Vector(const Vector &other)
    : Vector{other.size}
{
    for (auto i{0}; i < size; i+=4)
    {
        data[i] = other.data[i];
        data[i+1] = other.data[i+1];
        data[i+2] = other.data[i+2];
        data[i+3] = other.data[i+3];
    }
}

unsigned Vector::get_size() const
{
    return size;
}

double *Vector::get_data()
{
    return data;
}

double Vector::operator[](unsigned i) const
{
    return data[i];
}

double &Vector::operator[](unsigned i)
{
    return data[i];
}

double Vector::mean() const
{
    double sum{0};

    for (auto i{0}; i < size; i+=4)
    {
        sum += data[i];
        sum += data[i+1];
        sum += data[i+2];
        sum += data[i+3];
    }

    return sum / static_cast<double>(size);
}

double Vector::magnitude() 
{
    auto dot_prod{dot(*this)};
    return std::sqrt(dot_prod);
}

Vector& Vector::operator/(double div)
{

    __m256d divisor = _mm256_set1_pd(div);
    for (auto i{0}; i < size; i+=4)
    {
        __m256d values = _mm256_loadu_pd(&data[i]);
        __m256d result = _mm256_div_pd(values, divisor);
        _mm256_storeu_pd(&data[i],result);
    }

    return *this;
}

Vector& Vector::operator-(double sub)
{
    __m256d subtractor = _mm256_set1_pd(sub);

    for (auto i{0}; i < size; i+=4)
    {
        __m256d values = _mm256_loadu_pd(&data[i]);
        __m256d result = _mm256_sub_pd(values, subtractor);
        _mm256_storeu_pd(&data[i], result);
    }

    return *this;
}

double* Vector::getdata()
{
    return data;
}

double Vector::dot(Vector& rhs) const
{
    __m256d sum = _mm256_setzero_pd();

    for (auto i{0}; i < size; i+=4)
    {
        // Load 4 doubles from each vector into AVX registers
        __m256d vec1 = _mm256_loadu_pd(&data[i]);
        __m256d vec2 = _mm256_loadu_pd(&rhs.data[i]);

        // Perform element-wise multiplication
        __m256d prod = _mm256_mul_pd(vec1, vec2);

        // Accumulate the result into the sum
        sum = _mm256_add_pd(sum, prod);
    }
    double result[4];
    _mm256_storeu_pd(result, sum);
    double final_result = result[0] + result[1] + result[2] + result[3];

    return final_result;
}
