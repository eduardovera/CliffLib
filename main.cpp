#include <omp.h>
#include <iostream>
#include "multivector.hpp"
#include <dlib/image_io.h>
#include <dlib/external/libpng/png.h>

#include <chrono>
using namespace std;

using namespace CliffLib;

OrthonormalMetric<double> metric;

dlib::array2d<multivector<double>> double_to_multivec(const dlib::array2d<double> &m) {
    dlib::array2d<multivector<double>> m_;
    m_.set_size(m.nr(), m.nc());
    for (int j = 0, z = 1; j < m.nr(); j++) {
        for (int i = 0; i < m.nc(); i++, z++) {
            m_[j][i] = m[j][i] * e(z);
        }
    }
    return m_;
}

dlib::array2d<double> multivec_to_double(dlib::array2d<multivector<double>> &m) {
    dlib::array2d<double> m_;
    m_.set_size(m.nr(), m.nc());
    for (int j = 0; j < m.nr(); j++) {
        for (int i = 0; i < m.nc(); i++) {
            m_[j][i] = (m[j][i]).getScalar();
        }
    }
    return m_;
}


void flip_kernel(const dlib::array2d<double> &k, dlib::array2d<double> &k_) {
    for (int j = k_.nr() - 1, a = 0; j >= 0; j--, a++) {
        for (int i = k_.nc() - 1, b = 0; i >= 0; i--, b++) {
            k_[a][b] = k[j][i];
        }
    }
}

dlib::array2d<double> deconvolution(dlib::array2d<multivector<double>> &output, dlib::array2d<double> &K) {
    dlib::array2d<double> K_;
    K_.set_size(K.nr(), K.nc());
    flip_kernel(K, K_);

    dlib::array2d<multivector<double>> dims;
    dims.set_size(output.nr(), output.nc());
    for (int j = 0, z = 1; j < output.nr(); j++) {
        for (int i = 0; i < output.nc(); i++, z++) {
            dims[j][i] = e(z);
        }
    }

    int k_size = K_.nc() >> 1;
    dlib::array2d<double> input;
    input.set_size(output.nr(), output.nc());

    multivector<double> I;
    for (int j = 0; j < input.nr(); j = j+1) {
        for (int i = 0; i < input.nc(); i = i+1) {
            multivector<double> K_temp;
            for (int y = j - k_size, a = 0, it = 1; y <= j + k_size; y++, a++) {
                for (int x = i - k_size, b = 0; x <= i + k_size; x++, b++, it++) {
                    if (y < 0 || x < 0 || x >= output.nc() || y >= output.nr()) {
                        continue;
                    }
                    K_temp = K_temp + (K_[a][b] * dims[y][x]);
                    K_temp.handle_numeric_error();
                }
            }
            K_temp.handle_numeric_error();
            I = GP(IGP(output[j][i], K_temp, metric), dims[j][i], metric);
            I.handle_numeric_error();
            input[j][i] = I.getScalar();
        }
    }
    return input;
}

dlib::array2d<multivector<double>> convolution(dlib::array2d<multivector<double>> &I, dlib::array2d<double> &K) {
    dlib::array2d<double> K_;
    K_.set_size(K.nr(), K.nc());
    flip_kernel(K, K_);

    int k_size = K_.nc() >> 1;
    dlib::array2d<multivector<double>> output;
    output.set_size(I.nr(), I.nc());

    dlib::array2d<multivector<double>> dims;
    dims.set_size(I.nr(), I.nc());
    for (int j = 0, z = 1; j < I.nr(); j++) {
        for (int i = 0; i < I.nc(); i++, z++) {
            dims[j][i] = e(z);
        }
    }


    for (int j = 0; j < I.nr(); j++) {
        for (int i = 0; i < I.nc(); i++) {
            multivector<double> I_temp;
            multivector<double> K_temp;
            for (int y = j - k_size, a = 0; y <= j + k_size; y++, a++) {
                for (int x = i - k_size, b = 0; x <= i + k_size; x++, b++) {
                    if (y < 0 || x < 0 || x >= I.nc() || y >= I.nr()) {
                        continue;
                    }
                    I_temp = I_temp + I[y][x];
                    K_temp = K_temp + (K_[a][b] * dims[y][x]);
                }
            }
            I_temp.handle_numeric_error();
            K_temp.handle_numeric_error();
//            cout << I_temp << endl;
//            cout << K_temp << endl;
            output[j][i] = GP(I_temp, K_temp, metric);
            output[j][i].handle_numeric_error();
        }
    }
    return output;
}

void output_default_convolution(dlib::array2d<multivector<double>> &matrix, string filename) {
    cout << "Saving output... " << endl;
    auto sa_start = std::chrono::high_resolution_clock::now();
    dlib::array2d<unsigned char> output;
    output.set_size(matrix.nr(), matrix.nc());

    OrthonormalMetric<double> m;

    for (int j = 0; j < matrix.nr(); j++) {
        for (int i = 0; i < matrix.nc(); i++) {
            double g = ((matrix[j][i]).getScalar());

            g = g < 0 ? 0 : g;
            g = g > 255 ? 255 : g;

            output[j][i] = (unsigned char)(g);
        }
    }

    dlib::save_png(output, filename);
    auto sa_end = std::chrono::high_resolution_clock::now();
    auto sa_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(sa_end - sa_start).count()/1000.0;
    cout << "Done! (" << sa_elapsed << " seconds)" << endl;
}

template<class T>
void print(const dlib::array2d<T> &m) {
    for (int j = 0; j < m.nr(); j++) {
        for (int i = 0; i < m.nc(); i++) {
            cout << m[j][i] << endl;
        }
    }
}

dlib::array2d<double> conv(dlib::array2d<double> input, dlib::array2d<double> kernel, dlib::array2d<multivector<double>> &dict) {
    for (int j = 0, z = 1; j < dict.nr(); j++) {
        for (int i = 0; i < dict.nc(); i++, z++) {
            dict[j][i] = e(z);
        }
    }


}

int main() {

    CliffLib::N_DIMS = 50;

    dlib::array2d<double> img;

//    dlib::load_png(img, "/home/eduardovera/Workspace/CliffLib/6.png");

    img.set_size(3, 4);
    img[0][0] = 1;
    img[0][1] = 2;
    img[0][2] = 3;
    img[0][3] = 4;
    img[1][0] = 5;
    img[1][1] = 6;
    img[1][2] = 7;
    img[1][3] = 8;
    img[2][0] = 9;
    img[2][1] = 10;
    img[2][2] = 11;
    img[2][3] = 12;

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

    k2[0][0] = -1;
    k2[0][1] = -2;
    k2[0][2] = -1;
    k2[1][0] = 0;
    k2[1][1] = 0;
    k2[1][2] = 0;
    k2[2][0] = 1;
    k2[2][1] = 2;
    k2[2][2] = 1;

    dlib::array2d<multivector<double>> IMG = double_to_multivec(img);

    dlib::array2d<multivector<double>> O1 = convolution(IMG, k1);
    dlib::array2d<multivector<double>> O2 = convolution(O1, k2);
    dlib::array2d<double> I1 = deconvolution(O2, k2);
    dlib::array2d<multivector<double>> i = double_to_multivec(I1);
    dlib::array2d<double> I2 = deconvolution(i, k1);

    print(I2);



//    dlib::array2d<multivector<double>> K2 = double_to_multivec(k2);
//    dlib::array2d<multivector<double>> K = convolution(k1, K2);

//    K = double_to_multivec(multivec_to_double(K));

//    dlib::array2d<multivector<double>> IMG = convolution(img, K);
//    output_default_convolution(IMG, "comK1K2.png");
//    output_default_convolution(K);

//    dlib::array2d<multivector<double>> output = convolution(img, K);
//    dlib::array2d<double> k = multivec_to_double(K);

//    dlib::array2d<multivector<double>> b = deconvolution(IMG, K2);
//    dlib::array2d<multivector<double>> B = double_to_multivec(b);

//    dlib::array2d<multivector<double>> b = deconvolution(IMG, K2);

//    output_default_convolution(IMG, "semK2.png");

//    dlib::array2d<multivector<double>> t = deconvolution(b, K1);

}
