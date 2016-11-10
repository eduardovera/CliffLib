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

    OrthonormalMetric<double> m;
    multivector<double> k = UNDUAL(e(2), 3, m);
    cout << k << endl;
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

