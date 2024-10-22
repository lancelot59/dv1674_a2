/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "filters.hpp"
#include "matrix.hpp"
#include "ppm.hpp"
#include <cmath>
#include <pthread.h>


namespace Filter
{

    namespace Gauss
    {
        void get_weights(int n, double *weights_out)
        {
            for (auto i{0}; i <= n; i++)
            {
                double x{static_cast<double>(i) * max_x / n};
                weights_out[i] = exp(-x * x * pi);
            }
        }
    }

    // 6 Data needed by a thread to perform blur operation
    struct ThreadData {
        Matrix* matrix;
        unsigned char* R;
        unsigned char* G;
        unsigned char* B;
        double* w;
        int radius;
        unsigned x_size;
        unsigned y_size;
        int start;
        int end;
        ThreadData():
            matrix(nullptr),
            R(0),
            G(0),
            B(0),
            w(nullptr),
            radius(0),
            x_size(0),
            y_size(0),
            start(0),
            end(0) {}
        ThreadData(Matrix* matrix, unsigned char* R, unsigned char* G, unsigned char* B, double* weights, int radius,
            unsigned x_size, unsigned y_size, int start, int end):
        matrix(matrix), R(R), G(G), B(B), w(weights), radius(radius), x_size(x_size), y_size(y_size), start(start),
        end(end) {}
    };

    //Horizontal thread - operations on scratch
    void* horizontal_thread(void* arg) {
        auto* data = static_cast<ThreadData*>(arg);
        Matrix* matrix = data->matrix;
        unsigned char* R = data->R;
        unsigned char* G = data->G;
        unsigned char* B = data->B;
        double* w = data->w;
        const int radius = data->radius;
        const int x_size = data->x_size;
        const int y_size = data->y_size;
        const int start = data->start;
        const int end = data->end;

        //3 Moved getting first weight out of the loop - it does not change (read-only), and it is used to calculate
        // new values for colors of pixels
        auto w0 {w[0]};

        for (auto x{start}; x < end; x++)
        {
            for (auto y{0}; y < y_size; y++)
            {
                //4 Index is calculated before accessing the values in the arrays
                auto index {x + y * x_size};
                auto r{w0 * R[index]}, g{w0 * G[index]}, b{w0 * B[index]}, n{w0};

                for (auto wi{1}; wi <= radius; wi++) {
                    auto wc{w[wi]};

                    auto x2{x - wi};
                    if (x2 >= 0)
                    {
                        auto left_index {x2 + y * x_size};
                        r += wc * R[left_index];
                        g += wc * G[left_index];
                        b += wc * B[left_index];
                        n += wc;
                    }
                    x2 = x + wi;
                    if (x2 < x_size)
                    {
                        auto right_index {x2 + y * x_size};
                        r += wc * R[right_index];
                        g += wc * G[right_index];
                        b += wc * B[right_index];
                        n += wc;
                    }
                }
                //6 No mutex needed because each thread will access and change different pixel
                matrix->r(index) = r / n;
                matrix->g(index) = g / n;
                matrix->b(index) = b / n;

            }
        }

        return nullptr;
    }

    //Vertical thread - operations on dst
    void* vertical_thread(void* arg) {
        auto* data = static_cast<ThreadData*>(arg);
        Matrix* matrix = data->matrix;
        unsigned char* R = data->R;
        unsigned char* G = data->G;
        unsigned char* B = data->B;
        double* w = data->w;
        const int radius = data->radius;
        const int x_size = data->x_size;
        const int y_size = data->y_size;
        const int start = data->start;
        const int end = data->end;

        //3 Moved getting first weight out of the loop - it does not change (read-only), and it is used to calculate
        // new values for colors of pixels
        auto w0 {w[0]};

        for (auto y{start}; y < end; y++)
        {
            for (auto x{0}; x < x_size; x++)
            {
                //4 Index is calculated before accessing the values in the arrays
                auto index = x + y * x_size;
                if (index < 0 || index >= x_size * y_size) continue;
                auto r{w0 * R[index]}, g{w0 * G[index]}, b{w0 * B[index]}, n{w0};

                for (auto wi{1}; wi <= radius; wi++)
                {
                    auto wc{w[wi]};
                    auto y2{y - wi};

                    if (y2 >= 0)
                    {
                        auto left_index {x + y2 * x_size};
                        r += wc * R[left_index];
                        g += wc * G[left_index];
                        b += wc * B[left_index];
                        n += wc;
                    }
                    y2 = y + wi;
                    if (y2 < y_size)
                    {
                        auto right_index {x + y2 * x_size};
                        r += wc * R[right_index];
                        g += wc * G[right_index];
                        b += wc * B[right_index];
                        n += wc;
                    }
                }
                //6 No mutex needed because each thread will access and change different pixel
                matrix->r(index) = r / n;
                matrix->g(index) = g / n;
                matrix->b(index) = b / n;
            }
        }

        return nullptr;

    }

    //The blur function itself
    Matrix blur(Matrix m, const int radius, int num_threads)
    {
        Matrix scratch{PPM::max_dimension};
        auto dst{m};
        //1 Getting heights and width of the images before the loop, and then using it in conditions of for-loops
        auto x_size {dst.get_x_size()};
        auto y_size {dst.get_y_size()};
        //2 Gauss weight calculation is done before the loop - they are also read-only values used for calculations
        // of pixel's new colors values
        double w[Gauss::max_radius]{};
        Gauss::get_weights(radius, w);

        //5 Directly accessing R, G, and B - they are also read-only values used for calculations and this recduced
        // the amount of function calls
        auto dst_R {dst.get_R()}, dst_B {dst.get_B()}, dst_G {dst.get_G()};
        auto scr_R {scratch.get_R()}, scr_B {scratch.get_B()}, scr_G {scratch.get_G()};

        //6 Introducing multithreading
        pthread_t threads[num_threads];
        ThreadData thread_data[num_threads];
        int rows = x_size / num_threads;

        //Horizontal blurring
        for(int i = 0; i < num_threads; i++) {
            thread_data[i] = {
                &scratch,dst_R,dst_G,dst_B,w,radius,x_size,y_size,i * rows,
                static_cast<int>(i == num_threads - 1 ? x_size : (i + 1) * rows) //ensures all rows will be blurred
            };
            pthread_create(&threads[i], nullptr, horizontal_thread, &thread_data[i]);
        }

        for (int i = 0; i < num_threads; i++)
        {
            pthread_join(threads[i], nullptr);
        }

        //Vertical blurring
        for(int i = 0; i < num_threads; i++) {
            thread_data[i] = {
                &dst,scr_R,scr_G,scr_B,w,radius,x_size,y_size,i * rows,
                static_cast<int>(i == num_threads - 1 ? y_size : (i + 1) * rows) //ensures all rows will be blurred
            };
            pthread_create(&threads[i], nullptr, vertical_thread, &thread_data[i]);
        }

        for (int i = 0; i < num_threads; i++)
        {
            pthread_join(threads[i], nullptr);
        }


        return dst;
    }

}