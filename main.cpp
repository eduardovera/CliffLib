#include <omp.h>
#include <iostream>
#include "multivector.hpp"
#include <dlib/image_io.h>
#include <dlib/external/libpng/png.h>

#include <chrono>
using namespace std;

using namespace CliffLib;

dlib::array2d<multivector<double>> convolution(const dlib::array2d<unsigned char> &I, const dlib::array2d<double> &K) {
    OrthonormalMetric<double> m;
    dlib::array2d<multivector<double>> dims;
    dims.set_size(I.nr(), I.nc());
    for (int j = 0, z = 1; j < I.nr(); j++) {
        for (int i = 0; i < I.nc(); i++, z++) {
            dims[i][j] = MASKS[z];
//            cout << dims[i][j] << endl;
        }
    }

    dlib::array2d<multivector<double>> output;
    output.set_size(I.nr(), I.nc());

    int k_size = K.nc() >> 1;

    double scale = 1.0/255.0;

    int z = 0;
    #pragma omp parallel for
    for (int j = 0; j < I.nr(); j++) {
        for (int i = 0; i < I.nc(); i++) {
            multivector<double> I_temp;
            multivector<double> K_temp;
//            std::cout << "I: " << I_temp << std::endl;
//            std::cout << std::endl << "K: " << K_temp << std::endl;
            for (int y = j - k_size, a = 0; y <= j + k_size; y++, a++) {
                for (int x = i - k_size, b = 0; x <= i + k_size; x++, b++) {
//                    cout << "x: " << x << endl;
//                    cout << "y: " << y << endl;
                    if (y < 0 || x < 0 || x >= I.nc() || y >= I.nr()) {
                        continue;
                    }
//                    cout << "I(x, y): " << I[x][y] << endl;
//                    cout << "dims(x, y): " << dims[x][y] << endl;
//                    cout << "K(a, b): " << K[a][b] << endl;
//                    cout << "a: " << a << endl;
//                    cout << "b: " << b << endl;
                    I_temp = I_temp + (scale * (I[x][y]) * dims[x][y]);
                    K_temp = K_temp + (K[a][b] * dims[x][y]);
                }
            }
//            I_temp.handle_numeric_error();
//            K_temp.handle_numeric_error();
//            std::cout << "I: " << I_temp << std::endl;
//            std::cout << std::endl << "K: " << K_temp << std::endl;
            output[i][j] = GP(I_temp, K_temp, m);
            output[i][j].handle_numeric_error();
//            cout << "Computed " << ++z << " elements" << endl;
        }
    }
    return output;
}


int main() {

    CliffLib::N_DIMS = 400;

    CliffLib::build_masks();
    cout << "Loading image... " << endl;
    auto s_start = std::chrono::high_resolution_clock::now();
//    CliffLib::build_lookup_table(OrthonormalMetric<double>());

    dlib::array2d<dlib::rgb_pixel> imgRGB;
    dlib::load_image(imgRGB, "/home/eduardovera/Workspace/CliffLib/6.png");

    dlib::array2d<unsigned char> img;
    dlib::assign_image(img, imgRGB);

    dlib::array2d<double> kernel;

    kernel.set_size(3, 3);

    kernel[0][0] = 0;
    kernel[0][1] = -1;
    kernel[0][2] = 0;
    kernel[1][0] = -1;
    kernel[1][1] = 4;
    kernel[1][2] = -1;
    kernel[2][0] = 0;
    kernel[2][1] = -1;
    kernel[2][2] = 0;

    auto s_end = std::chrono::high_resolution_clock::now();
    auto s_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(s_end - s_start).count()/1000.0;
    cout << "Done! (" << s_elapsed << " seconds)" << endl;

    cout << "Starting convolution... " << endl;
    auto c_start = std::chrono::high_resolution_clock::now();
    dlib::array2d<multivector<double>> G = convolution(img, kernel);
    auto c_end = std::chrono::high_resolution_clock::now();
    auto c_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(c_end - c_start).count()/1000.0;
    cout << "Done! (" << c_elapsed << " seconds)" << endl;

    cout << "Saving output... " << endl;
    auto sa_start = std::chrono::high_resolution_clock::now();
    dlib::array2d<unsigned char> output;
    output.set_size(img.nr(), img.nc());

    OrthonormalMetric<double> m;

    for (int j = 0; j < img.nr(); j++) {
        for (int i = 0; i < img.nc(); i++) {
            double g = (G[i][j]).getCoeff();
            output[i][j] = (unsigned char)(255 * g);
        }
    }

    dlib::save_png(output, "output.png");
    auto sa_end = std::chrono::high_resolution_clock::now();
    auto sa_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(sa_end - sa_start).count()/1000.0;
    cout << "Done! (" << sa_elapsed << " seconds)" << endl;

}

