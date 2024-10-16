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
    for (int i = 0; i < size; i+=4)
    {
        __m256d dataVec = _mm256_loadu_pd(&other.data[i]);
        _mm256_storeu_pd(&data[i], dataVec);
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
    __m256d sumVec = _mm256_setzero_pd();
    int i = 0;

    for (; i <= size -4; i+=4)
    {
        __m256d dataVec = _mm256_loadu_pd(&data[i]);
        sumVec = _mm256_add_pd(sumVec, dataVec);
    }
    double sum = 0;
    double temp[4];
    _mm256_storeu_pd(temp, sumVec);
    sum += temp[0] + temp[1] + temp[2] + temp[3];
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

double Vector::dot(Vector& rhs) const
{
    __m256d sum = _mm256_setzero_pd();

    for(int i = 0; i < size; i+=4)
    {
        __m256d a = _mm256_loadu_pd(&data[i]);
        __m256d b = _mm256_loadu_pd(&rhs[i]);

        __m256d product = _mm256_mul_pd(a,b);

        sum = _mm256_add_pd(sum,product);
    }

    double result[4];
    _mm256_storeu_pd(result,sum);
    return result[0] + result[1] + result[2] + result[3];
}