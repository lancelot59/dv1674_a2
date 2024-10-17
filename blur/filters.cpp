/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "filters.hpp"
#include "matrix.hpp"
#include "ppm.hpp"
#include <cmath>
#include <pthread.h>
#include <vector>


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

    struct ThreadData {
        Matrix* dst;
        Matrix* scratch;
        double* w;
        int radius;
        unsigned x_size;
        unsigned y_size;
        int start;
        int end;
        ThreadData()
            : dst(nullptr), scratch(nullptr), w(nullptr), radius(0),
            x_size(0), y_size(0), start(0), end(0) {}
        ThreadData(Matrix* d, Matrix* s, double* weights, int r, unsigned xs, unsigned ys, int st, int en)
            : dst(d), scratch(s), w(weights), radius(r), x_size(xs), y_size(ys), start(st), end(en) {}
    };

    void* thread1(void* arg) {
        ThreadData* data = static_cast<ThreadData*>(arg);
        Matrix* dst = data->dst;
        Matrix* scratch = data->scratch;
        double* w = data->w;
        const int radius = data->radius;
        unsigned x_size = data->x_size;
        unsigned y_size = data->y_size;
        int start = data->start;
        int end = data->end;

        // 3 move this outside of the loop
        auto w0 {w[0]};

        //4 Moved those calling for those arrays outside of the loop. They are called once instead of multiple times
        auto R_dst = dst->get_R(), G_dst = dst->get_G(), B_dst = dst->get_B();
        auto R_scratch = scratch->get_R(), G_scratch = scratch->get_G(), B_scratch = scratch->get_B();


         for (auto x{start}; x < end; x++)
        {
            for (auto y{0}; y < y_size; y++)
            {
                auto index {y * x_size + x};
                auto r{w0 * R_dst[index]}, g{w0 * G_dst[index]}, b{w0 * B_dst[index]}, n{w0};

                for (auto wi{1}; wi <= radius; wi++)
                {
                    auto wc{w[wi]};
                    auto x2{x - wi};

                    if (x2 >= 0)
                    {
                        auto left_index = y * x_size + x2;
                        r += wc * R_dst[left_index];
                        g += wc * G_dst[left_index];
                        b += wc * B_dst[left_index];
                        n += wc;
                    }
                    x2 = x + wi;
                    if (x2 < x_size)
                    {
                        auto right_index = y * x_size + x2;
                        r += wc * R_dst[right_index];
                        g += wc * G_dst[right_index];
                        b += wc * B_dst[right_index];
                        n += wc;
                    }
                }
                // those one stays because R, G, B are const
                R_scratch[index] = r / n;
                G_scratch[index] = g / n;
                B_scratch[index] = b / n;
            }
        }

        for (auto x{start}; x < end; x++)
        {
            for (auto y{0}; y < y_size; y++)
            {
                auto index {y * x_size + x};
                auto r{w0 * R_scratch[index]}, g{w0 * G_scratch[index]}, b{w0 * B_scratch[index]}, n{w0};

                for (auto wi{1}; wi <= radius; wi++)
                {
                    auto wc{w[wi]};
                    auto y2{y - wi};

                    if (y2 >= 0)
                    {
                        auto left_index = y2 * x_size + x;
                        r += wc * R_scratch[left_index];
                        g += wc * G_scratch[left_index];
                        b += wc * B_scratch[left_index];
                        n += wc;
                    }
                    y2 = y + wi;
                    if (y2 < y_size)
                    {
                        auto right_index = y2 * x_size + x;
                        r += wc * R_scratch[right_index];
                        g += wc * G_scratch[right_index];
                        b += wc * B_scratch[right_index];
                        n += wc;
                    }
                }
                R_dst[index] = r / n;
                G_dst[index] = g / n;
                B_dst[index] = b / n;
            }
        }
        return nullptr;

    }

    Matrix blur(Matrix m, const int radius)
    {
        Matrix scratch{PPM::max_dimension};
        auto dst{m};

        // 1 adding getting date before the loop
        auto x_size {dst.get_x_size()};
        auto y_size {dst.get_y_size()};

        // 2 Gauss operation before the loop not inside
        double w[Gauss::max_radius]{};
        Gauss::get_weights(radius, w);

        //5 removed the const from getter for R, G, and B in matrix.cpp
        const int num_threads = 4;
        pthread_t threads[num_threads];
        std::vector<ThreadData> thread_data(num_threads);

        int rows = x_size / num_threads;

        for(int i = 0; i < num_threads; i++) {
            thread_data[i] = {&scratch, &dst, w, radius, x_size, y_size, i * rows, (i + 1) * rows};
            if (i == num_threads - 1) {
                thread_data[i].end = y_size;
            }
            pthread_create(&threads[i], nullptr, thread1, &thread_data[i]);
        }

        for (int i = 0; i < num_threads; i++) {
            pthread_join(threads[i], nullptr);
        }

        return dst;
    }

}
