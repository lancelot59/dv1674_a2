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

void Vector::index(unsigned i, double* returnarray)
{
    returnarray[0] = data[i];
    returnarray[1] = data[i+1];
    returnarray[2] = data[i+2];
    returnarray[3] = data[i+3];
}

double Vector::dot(Vector& rhs) const
{
    double result{0};
    double array[4];


    for (auto i{0}; i < size; i+=4)
    {
        rhs.index(i, array);
        result += data[i] * array[0];
        result += data[i+1] * array[1];
        result += data[i+2] * array[2];
        result += data[i+3] * array[3];
    }

    return result;
}