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


void flip_kernel(const dlib::array2d<multivector<double>> &k, dlib::array2d<multivector<double>> &k_) {
    for (int j = k_.nr() - 1, a = 0; j >= 0; j--, a++) {
        for (int i = k_.nc() - 1, b = 0; i >= 0; i--, b++) {
            k_[a][b] = k[j][i];
        }
    }
}

dlib::array2d<multivector<double>> deconvolution(dlib::array2d<multivector<double>> &output, dlib::array2d<multivector<double>> &K) {
    dlib::array2d<multivector<double>> K_;
    K_.set_size(K.nr(), K.nc());
    flip_kernel(K, K_);


    multivector<double> K_temp;

    for (int j = 0; j < K_.nr(); j++) {
        for (int i = 0; i < K_.nc(); i++) {
            K_temp = K_temp + K_[j][i];
        }
    }
    K_temp.handle_numeric_error();

    int k_size = K_.nc() >> 1;
    dlib::array2d<multivector<double>> input;
    input.set_size(output.nr() - (2 * k_size), output.nc() - (2 * k_size));

    multivector<double> I;

    for (int j = k_size, r = 0; j < output.nr() - k_size; j = j+1, r++) {
        for (int i = k_size, s = 0; i < output.nc() - k_size; i = i+1, s++) {

            if ((j >= output.nr() - (2 * k_size)) || (i >= output.nc() - (2 * k_size))) {
//                continue;
            }
            I = IGP(output[j][i], K_temp, metric);
            I.handle_numeric_error();

//            for (int y = 0, it = 1; y < input.nr(); y++) {
//                for (int x = 0; x < input.nc(); x++, it++) {
//                    input[r][s] = I.get_portion(it);
//                }
//            }

            cout << I << endl;

            getchar();

//            for (int y = j - k_size, a = 0, it = 1; y <= j + k_size; y++, a++) {
//                for (int x = i - k_size, b = 0; x <= i + k_size; x++, b++, it++) {
//                    if (x >= input.nc() || y >= input.nr() || x < 0 || y < 0) {
//                        continue;
//                    }
//                    input[y][x] = I.get_portion(it);
//                }
//            }

        }
    }

    for (int j = 0; j < input.nr(); j++) {
        for (int i = 0; i < input.nc(); i++) {
            cout << input[j][i] << endl;
        }
    }




}

dlib::array2d<multivector<double>> convolution(const dlib::array2d<double> &I, const dlib::array2d<multivector<double>> &K) {
    dlib::array2d<multivector<double>> K_;
    K_.set_size(K.nr(), K.nc());
    flip_kernel(K, K_);

    int k_size = K_.nc() >> 1;
    dlib::array2d<multivector<double>> output;
    output.set_size(I.nr() + (2 * k_size), I.nc() + (2 * k_size));

    dlib::array2d<multivector<double>> dims;
    dims.set_size(K.nr(), K.nc());
    for (int j = 0, z = 1; j < K.nr(); j++) {
        for (int i = 0; i < K.nc(); i++, z++) {
            dims[j][i] = e(z);
        }
    }

    for (int j = 0, im_j = 0; j < output.nr(); j++, im_j++) {
        for (int i = 0, im_i = 0; i < output.nc(); i++, im_i++) {
            multivector<double> I_temp;
            multivector<double> K_temp;
            for (int y = j - k_size, a = 0; y <= j + k_size; y++, a++) {
                for (int x = i - k_size, b = 0; x <= i + k_size; x++, b++) {
                    if (y < 0 || x < 0 || x >= output.nc() || y >= output.nr()) {
                        continue;
                    }
                    double img_at;
                    if (y - k_size < 0 || x - k_size < 0 || y - k_size + 1 > I.nr() || x - k_size + 1> I.nc()) {
                        img_at = 0;
                    } else {
                        img_at = I[y - k_size][x - k_size];
                    }
                    I_temp = I_temp + (img_at * dims[a][b]);
                    K_temp = K_temp + (K_[a][b]);
                }
            }
            I_temp.handle_numeric_error();
            K_temp.handle_numeric_error();
            output[j][i] = GP(I_temp, K_temp, metric);
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
            double g = ((matrix[j][i]).getScalar());

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

    CliffLib::N_DIMS = 50;

//    ConformalMetric<double> mC;

//    multivector<double> A = e(0)^e(20)^e(50);
////    multivector<double> B = GP(GP(A, e(4)^e(50), mC), INV(A, mC), mC);
//    cout << GP(A, e(0)^e(50), mC) << endl;
//    return 0;

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

    k2[0][0] = 1;
    k2[0][1] = 2;
    k2[0][2] = 1;
    k2[1][0] = 0;
    k2[1][1] = 0;
    k2[1][2] = 0;
    k2[2][0] = -1;
    k2[2][1] = -2;
    k2[2][2] = -1;

    dlib::array2d<multivector<double>> K1 = double_to_multivec(k1);
    dlib::array2d<multivector<double>> K2 = double_to_multivec(k2);

    dlib::array2d<multivector<double>> K = convolution(img, K2);
    dlib::array2d<multivector<double>> K_INV = deconvolution(K, K2);


//    dlib::array2d<multivector<double>> IMG = convolution(img, K);
//    output_default_convolution(K);

}

