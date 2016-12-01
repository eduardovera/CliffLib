#include <omp.h>
#include <iostream>
#include "multivector.hpp"
#include <dlib/image_io.h>
#include <dlib/external/libpng/png.h>

#include <chrono>
using namespace std;

using namespace CliffLib;

void flip_kernel(const dlib::array2d<double> &k, dlib::array2d<double> &k_) {
//    k_.set_size(k.nr(), k.nc());
    for (int j = k_.nr() - 1, a = 0; j >= 0; j--, a++) {
        for (int i = k_.nc() - 1, b = 0; i >= 0; i--, b++) {
            k_[a][b] = k[j][i];
        }
    }
}

dlib::array2d<multivector<double>> convolution(const dlib::array2d<double> &I, dlib::array2d<double> &K) {
    OrthonormalMetric<double> m;
    dlib::array2d<multivector<double>> dims;
    dims.set_size(I.nr(), I.nc());
    for (int j = 0, z = 1; j < K.nr(); j++) {
        for (int i = 0; i < K.nc(); i++, z++) {
            dims[j][i] = MASKS[z];
        }
    }

    dlib::array2d<double> K_;
    K_.set_size(K.nr(), K.nc());
    flip_kernel(K, K_);

    dlib::array2d<multivector<double>> output;
    output.set_size(I.nr(), I.nc());

    int k_size = K_.nc() >> 1;

    for (int j = 0; j < I.nr(); j++) {
        for (int i = 0; i < I.nc(); i++) {
            multivector<double> I_temp;
            multivector<double> K_temp;
            for (int y = j - k_size, a = 0; y <= j + k_size; y++, a++) {
                for (int x = i - k_size, b = 0; x <= i + k_size; x++, b++) {
                    if (y < 0 || x < 0 || x >= I.nc() || y >= I.nr()) {
                        continue;
                    }
                    I_temp = I_temp + (I[y][x] * dims[a][b]);
                    K_temp = K_temp + (K_[a][b] * dims[a][b]);
                }
            }
            I_temp.handle_numeric_error();
            K_temp.handle_numeric_error();
            output[j][i] = GP(I_temp, K_temp, m);
            output[j][i].handle_numeric_error();
//            cout << "I: " << I_temp << endl;
//            cout << "K: " << K_temp << endl;
//            cout << "OUT: " << output[j][i] << endl;
//            getchar();
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

    CliffLib::N_DIMS = 9;

    CliffLib::build_masks();
    cout << "Loading image... " << endl;
    auto s_start = std::chrono::high_resolution_clock::now();
//    CliffLib::build_lookup_table(OrthonormalMetric<double>());

    dlib::array2d<double> img;
    img.set_size(3, 3);

    img[0][0] = 1;
    img[0][1] = 2;
    img[0][2] = 3;
    img[1][0] = 4;
    img[1][1] = 5;
    img[1][2] = 6;
    img[2][0] = 7;
    img[2][1] = 8;
    img[2][2] = 9;


//    dlib::load_png(img, "/home/eduardovera/Workspace/CliffLib/6.png");
    dlib::array2d<double> k1;
    k1.set_size(3, 3);

    k1[0][0] = -1;
    k1[0][1] = -1;
    k1[0][2] = -1;
    k1[1][0] = -1;
    k1[1][1] = 8;
    k1[1][2] = -1;
    k1[2][0] = -1;
    k1[2][1] = -1;
    k1[2][2] = -1;


    dlib::array2d<double> k2;
    k2.set_size(3, 3);

    k2[0][0] = 1;
    k2[0][1] = 2;
    k2[0][2] = 1;
    k2[1][0] = 0;
    k2[1][1] = 0;
    k2[1][2] = 0;
    k2[2][0] = -1;
    k2[2][1] = -2;
    k2[2][2] = -1;


    auto s_end = std::chrono::high_resolution_clock::now();
    auto s_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(s_end - s_start).count()/1000.0;
    cout << "Done! (" << s_elapsed << " seconds)" << endl;

    cout << "Starting convolution... " << endl;
    auto c_start = std::chrono::high_resolution_clock::now();
    dlib::array2d<multivector<double>> temp = convolution(img, k1);
    auto c_end = std::chrono::high_resolution_clock::now();
    auto c_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(c_end - c_start).count()/1000.0;
    cout << "Done! (" << c_elapsed << " seconds)" << endl;

    dlib::array2d<double> I_;
    I_.set_size(temp.nr(), temp.nc());

    for (int j = 0; j < temp.nr(); j++) {
        for (int i = 0; i < temp.nc(); i++) {
//            cout << (temp[j][i]).getCoeff() << endl;
            I_[j][i] = (temp[j][i]).getCoeff();
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



//    G = convolution(img, I_);

//    dlib::array2d<double> temp;
//    temp.set_size(G.nr(), G.nc());

    dlib::array2d<multivector<double>> G = convolution(I_, k2);
    for (int j = 0; j < temp.nr(); j++) {
        for (int i = 0; i < temp.nc(); i++) {
            cout << (G[j][i]).getCoeff() << endl;
        }
    }

    output_default_convolution(G);


}

