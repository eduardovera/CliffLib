#include <iostream>
#include "multivector.hpp"
#include <dlib/image_io.h>
#include <dlib/external/libpng/png.h>

using namespace std;

using namespace CliffLib;

dlib::array2d<multivector<double>> convolution(const dlib::array2d<double> &I, const dlib::array2d<double> &K) {
    OrthonormalMetric<double> m;
    dlib::array2d<multivector<double>> dims;
    dims.set_size(I.nr(), I.nc());
    for (int j = 0, z = 1; j < I.nr(); j++) {
        for (int i = 0; i < I.nc(); i++, z++) {
            dims[i][j] = e(z);
//            cout << dims[i][j] << endl;
        }
    }

    dlib::array2d<multivector<double>> output;
    output.set_size(I.nr(), I.nc());

    int k_size = K.nc() >> 1;

    for (int j = 0; j < I.nr(); j++) {
        for (int i = 0; i < I.nc(); i++) {
            multivector<double> I_temp = SCALAR<double>();
            multivector<double> K_temp = SCALAR<double>();
            std::cout << "I: " << I_temp << std::endl;
            std::cout << std::endl << "K: " << K_temp << std::endl;
            for (int y = j - k_size, a = 0; y <= j + k_size; y++, a++) {
                for (int x = i - k_size, b = 0; x <= i + k_size; x++, b++) {
//                    cout << "x: " << x << endl;
//                    cout << "y: " << y << endl;
                    if (y < 0 || x < 0 || x >= I.nc() || y >= I.nr()) {
                        continue;
                    }
                    cout << "I(x, y): " << I[x][y] << endl;
                    cout << "dims(x, y): " << dims[x][y] << endl;
                    cout << "K(a, b): " << K[a][b] << endl;
                    cout << "a: " << a << endl;
                    cout << "b: " << b << endl;
                    I_temp = I_temp + (I[x][y] * dims[x][y]);
                    K_temp = K_temp + (K[x - i + k_size][y - j + k_size] * dims[x][y]);
                    I_temp.handle_numeric_error();
                    K_temp.handle_numeric_error();
                }
            }
            std::cout << "I: " << I_temp << std::endl;
            std::cout << std::endl << "K: " << K_temp << std::endl;
            output[i][j] = GP(I_temp, K_temp, m);
            output[i][j].handle_numeric_error();
        }
    }
    return output;
}


int main() {

    CliffLib::N_DIMS = 10000;

    dlib::array2d<double> img;
    dlib::load_png(img, "/home/eduardovera/Workspace/CliffLib/6.png");

    dlib::array2d<double> prewitt_operator_X;
    dlib::array2d<double> prewitt_operator_Y;

    prewitt_operator_X.set_size(3, 3);
    prewitt_operator_Y.set_size(3, 3);

    prewitt_operator_X[0][0] = -1;
    prewitt_operator_X[1][0] = -1;
    prewitt_operator_X[2][0] = -1;
    prewitt_operator_X[0][2] = +1;
    prewitt_operator_X[1][2] = +1;
    prewitt_operator_X[2][2] = +1;

    prewitt_operator_Y[0][0] = -1;
    prewitt_operator_Y[0][1] = -1;
    prewitt_operator_Y[0][2] = -1;
    prewitt_operator_Y[2][0] = +1;
    prewitt_operator_Y[2][1] = +1;
    prewitt_operator_Y[2][2] = +1;

    dlib::array2d<multivector<double>> Gx = convolution(img, prewitt_operator_X);
//    dlib::array2d<multivector<double>> Gy = convolution(img, prewitt_operator_Y);

}

