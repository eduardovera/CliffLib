#include <omp.h>
#include <iostream>
#include "multivector.hpp"
#include <dlib/image_io.h>
#include <dlib/external/libpng/png.h>

#include <chrono>
using namespace std;

using namespace CliffLib;

dlib::array2d<multivector<double>> convolution(const dlib::array2d<double> &I, const dlib::array2d<double> &K) {
    OrthonormalMetric<double> m;
    dlib::array2d<multivector<double>> dims;
    dims.set_size(I.nr(), I.nc());
    for (int j = 0, z = 1; j < I.nr(); j++) {
        for (int i = 0; i < I.nc(); i++, z++) {
            dims[j][i] = MASKS[z];
        }
    }

    dlib::array2d<multivector<double>> output;
    output.set_size(I.nr(), I.nc());

    int k_size = K.nc() >> 1;

    for (int j = 0; j < I.nr(); j++) {
        for (int i = 0; i < I.nc(); i++) {
            multivector<double> I_temp;
            multivector<double> K_temp;
            for (int y = j - k_size, a = 0; y <= j + k_size; y++, a++) {
                for (int x = i - k_size, b = 0; x <= i + k_size; x++, b++) {
                    if (y < 0 || x < 0 || x >= I.nc() || y >= I.nr()) {
                        continue;
                    }
                    I_temp = I_temp + (((I[y][x])) * dims[y][x]);
                    K_temp = K_temp + (K[a][b] * dims[y][x]);
//                    cout << "I: " << I[y][x] << endl;
//                    cout << "K: " << K[a][b] << endl;
//                    getchar();
                }
            }
            I_temp.handle_numeric_error();
            K_temp.handle_numeric_error();
            output[j][i] = GP(I_temp, K_temp, m);
            output[j][i].handle_numeric_error();
        }
    }
    return output;
}

void output_default_convolution(dlib::array2d<multivector<double>> &matrix) {
    cout << "Saving output... " << endl;
    auto sa_start = std::chrono::high_resolution_clock::now();
    dlib::array2d<unsigned char> output;
    output.set_size(matrix.nr(), matrix.nc());

    OrthonormalMetric<double> m;

    for (int j = 0; j < matrix.nr(); j++) {
        for (int i = 0; i < matrix.nc(); i++) {
            double g = ((matrix[j][i]).getCoeff());

            g = g < 0 ? 0 : g;
            g = g > 255 ? 255 : g;

            output[j][i] = (unsigned char)(g);
        }
    }

    dlib::save_png(output, "output.png");
    auto sa_end = std::chrono::high_resolution_clock::now();
    auto sa_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(sa_end - sa_start).count()/1000.0;
    cout << "Done! (" << sa_elapsed << " seconds)" << endl;
}


int main() {

    CliffLib::N_DIMS = 300;

    CliffLib::build_masks();
    cout << "Loading image... " << endl;
    auto s_start = std::chrono::high_resolution_clock::now();
//    CliffLib::build_lookup_table(OrthonormalMetric<double>());

    dlib::array2d<double> img;
    dlib::load_png(img, "/home/eduardovera/Workspace/CliffLib/3x3.png");
    dlib::array2d<double> k1;
    k1.set_size(3, 3);

    k1[0][0] = 1;
    k1[0][1] = 2;
    k1[0][2] = 3;
    k1[1][0] = 4;
    k1[1][1] = 5;
    k1[1][2] = 6;
    k1[2][0] = -38;
    k1[2][1] = 8;
    k1[2][2] = 9;


    dlib::array2d<double> k2;
    k2.set_size(3, 3);

    k2[0][0] = -1;
    k2[0][1] = -1;
    k2[0][2] = -1;
    k2[1][0] = -1;
    k2[1][1] = 8;
    k2[1][2] = -1;
    k2[2][0] = -1;
    k2[2][1] = -1;
    k2[2][2] = -1;


    auto s_end = std::chrono::high_resolution_clock::now();
    auto s_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(s_end - s_start).count()/1000.0;
    cout << "Done! (" << s_elapsed << " seconds)" << endl;

    cout << "Starting convolution... " << endl;
    auto c_start = std::chrono::high_resolution_clock::now();
    dlib::array2d<multivector<double>> G = convolution(k1, k2);
    auto c_end = std::chrono::high_resolution_clock::now();
    auto c_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(c_end - c_start).count()/1000.0;
    cout << "Done! (" << c_elapsed << " seconds)" << endl;

    dlib::array2d<double> I_;
    I_.set_size(3, 3);

    for (int j = 0; j < G.nr(); j++) {
        for (int i = 0; i < G.nc(); i++) {
            I_[j][i] = (G[j][i]).getCoeff();
//            cout << (G[j][i]).getCoeff() << endl;
        }
    }

//    dlib::array2d<double> k;
//    k.set_size(3, 3);

//    k[0][0] = 4;
//    k[0][1] = 0;
//    k[0][2] = -4;
//    k[1][0] = 0;
//    k[1][1] = 0;
//    k[1][2] = 0;
//    k[2][0] = -4;
//    k[2][1] = 0;
//    k[2][2] = 4;



    G = convolution(img, I_);

//    dlib::array2d<double> temp;
//    temp.set_size(G.nr(), G.nc());

    for (int j = 0; j < G.nr(); j++) {
        for (int i = 0; i < G.nc(); i++) {
//            I_[j][i] = (G[j][i]).getCoeff();
            cout << (G[j][i]).getCoeff() << endl;
        }
    }

//    G = convolution(temp, k2);
//    output_default_convolution(G);


}

