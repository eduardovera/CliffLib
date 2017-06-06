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
            m_[j][i] = scalar(m[j][i]);
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

dlib::array2d<multivector<double>> deconvolution(dlib::array2d<multivector<double>> &output, dlib::array2d<double> &K, int seed) {
    dlib::array2d<double> K_;
    K_.set_size(K.nr(), K.nc());
    flip_kernel(K, K_);

    dlib::array2d<multivector<double>> dims;
    dims.set_size(output.nr(), output.nc());
    for (int j = 0, z = seed; j < output.nr(); j++) {
        for (int i = 0; i < output.nc(); i++, z++) {
            dims[j][i] = e(z);
        }
    }

    int k_size = K_.nc() >> 1;
    dlib::array2d<multivector<double>> input;
    input.set_size(output.nr(), output.nc());

    multivector<double> I;
    for (int j = 0; j < input.nr(); j = j+1) {
        for (int i = 0; i < input.nc(); i = i+1) {
            multivector<double> K_temp;
            for (int y = j - k_size, a = 0, it = 1; y <= j + k_size; y++, a++) {
                for (int x = i - k_size, b = 0; x <= i + k_size; x++, b++, it++) {
                    if (y < 0 || x < 0 || x >= input.nc() || y >= input.nr()) {
                        continue;
                    }
//                    cout << K_[a][b] << " " << dims[y][x] << endl;
//                    getchar();
                    K_temp = K_temp + (K_[a][b] * dims[y][x]);
                }
            }
            K_temp.handle_numeric_error();
//            cout << output[j][i] << endl;
//            cout << K_temp << endl;
//            cout << IGP(output[j][i], K_temp, metric) << endl;
//            getchar();
            I = IGP(IGP(output[j][i], K_temp, metric), dims[j][i], metric);
            I.handle_numeric_error();
            input[j][i] = I;
//            cout << K_temp << endl;
//            cout << IGP(IGP(output[j][i], K_temp, metric), dims[j][i], metric) << endl;
//            getchar();
//            cout << IGP(I, dims[j][i], metric).getScalar() << endl;
        }
    }
    return input;
}

dlib::array2d<multivector<double>> convolution(dlib::array2d<multivector<double>> &I, dlib::array2d<double> &K, int seed) {
    dlib::array2d<double> K_;
    K_.set_size(K.nr(), K.nc());
    flip_kernel(K, K_);

    int k_size = K_.nc() >> 1;
    dlib::array2d<multivector<double>> output;
    output.set_size(I.nr(), I.nc());

    dlib::array2d<multivector<double>> dims;
    dims.set_size(I.nr(), I.nc());
    for (int j = 0, z = seed; j < I.nr(); j++) {
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
                    I_temp = I_temp + GP(I[y][x], dims[y][x], metric);
                    K_temp = K_temp + (K_[a][b] * dims[y][x]);
                }
            }
            I_temp.handle_numeric_error();
            K_temp.handle_numeric_error();
//            cout << I_temp << endl;
//            cout << K_temp << endl;
//            cout << "CONV: " << K_temp << endl;
//            cout << "I (conv): " << I_temp << endl;
            output[j][i] = GP(I_temp, K_temp, metric);
            output[j][i].handle_numeric_error();
//            cout << output[j][i] << endl;
//            getchar();
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

    CliffLib::N_DIMS = 12;

//    build_lookup_table(metric);

    multivector<double> I = 0*e(1)+2*e(2)+3*e(3)+5*e(4)+0*e(5);
    multivector<double> I_ = 4*e(1)+16*e(2)+29*e(3)+31*e(4)+10*e(5);

    multivector<double> VV = IGP(I_, I, metric);
    VV.handle_numeric_error();

    multivector<double> B = (I^I_);
    double n_b = 1.0 / SQR_NORM_REVERSE(B, metric);
    B = B * n_b;

    double norm = -.210526/9.68382e-05;

    cout << (B) << endl;

    cout << (VV) << endl;

    double angle = 0.5 * atan2(norm, 0);

    cout << angle << endl;

    multivector<double> V = cos(angle) - (sin(angle) * B);
    V.handle_numeric_error();

    multivector<double> rebuild = IGP(GP(V, I, metric), V, metric);
    rebuild.handle_numeric_error();
    cout << rebuild << endl;

}
