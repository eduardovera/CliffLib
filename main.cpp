#include <iostream>
#include "multivector.hpp"
#include <dlib/image_io.h>
#include <dlib/external/libpng/png.h>

using namespace std;

using namespace CliffLib;


int main() {

    CliffLib::N_DIMS = 10000;

    dlib::array2d<double> img;
    dlib::load_png(img, "/home/eduardovera/Workspace/CliffLib/input.png");

    dlib::array2d<multivector<double>> blades;
    blades.set_size(img.nr(), img.nc());

    dlib::array2d<multivector<double>> input_img;
    input_img.set_size(img.nr(), img.nc());


    int k = 1;
    for (int j = 0; j < img.nr(); ++j) {
        for (int i = 0; i < img.nc(); ++i) {
            blades[i][j] = e(k);
            input_img[i][j] = img[i][j] * blades[i][j];
            ++k;
        }
    }

    OrthonormalMetric<double> m;
    cout << GP(blades[45][32], input_img[45][32], m) << endl;


//    multivector<double> I = (3*e(1)) + (5*e(2)) + (4*e(4)) + (9*e(5));
//    multivector<double> K1 = (5*e(1)) + (6*e(2)) + (8*e(4)) + (9*e(5));
//    multivector<double> K2 = (3*e(1)) + (5*e(2)) + (3*e(4)) + (1*e(5));

//    multivector<double> R = GP(GP(I, K1, m), K2, m);

//    multivector<double> K = GP(K1, K2, m);

//    multivector<double> I_1 = IGP(R, K, m);
//    multivector<double> I_2 = IGP(R, K2, m);
//    I_2 = IGP(I_2, K1, m);

//    cout << I_1 << endl;
//    cout << I_2 << endl;


}

