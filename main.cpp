#include <iostream>
#include "multivector.hpp"

using namespace std;

using namespace CliffLib;


int main() {

    /**
     *
     * TO DO:
     *
     * - Change pointers to direct access
     * - Make orthogonal metric accept (P,Q,R) pattern
     *
     **/

    CliffLib::N_DIMS = 5;

    OrthonormalMetric<double> m;
    multivector<double> k1 = 5 + 2*(e(2)^e(3)) + 4*(e(1)^e(2)) + 0.5*(e(3)^e(1));
    multivector<double> k2 = e(1)^e(2)^e(3);
    cout << k1 << endl;
    cout << k1.get_type(m) << endl;
//    double scaling;
//    std::vector<multivector<double>> k = FACTORIZE(((10*e(1)^e(2)^e(5))+(3*e(1)^e(3)^e(5))+(4*(e(5)^e(3)^e(2)))), m, scaling);
//    multivector<double> check = scaling * k[0];
//    for (int i = 0; i < k.size(); i++) {
//        cout << scaling * k[i] << endl;
//        if (i >= 1) {
//            check = check ^ (scaling * k[i]);
//        }
//    }

//    cout << scaling * check << endl;

}

