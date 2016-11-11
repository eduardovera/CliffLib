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

    CliffLib::N_DIMS = 4;

    OrthonormalMetric<double> m;
    multivector<double> multivec = GP(e(1)^e(2), e(2)^e(4), m);
//    double k;
//    for (auto f : FACTORIZE(multivec, m, k)) {
//        cout << f << endl;
//    }
    cout << multivec << endl;
    cout << multivec.get_type(m) << endl;

//    multivector<double> multivec = GP((e(1)^e(2))+(e(1)^e(3)), (e(1)^e(2)) +(e(1)^e(3)) , m);
//    cout << multivec << endl;

//    cout << k1 << endl;
//    cout << k1.get_type(m) << endl;
//    double scaling;
//    std::vector<multivector<double>> k = FACTORIZE(((3*e(1)^e(2))+(4*e(1)^e(3))), m, scaling);
//    cout << scaling << endl;
//    multivector<double> check = SCALAR<double>();
//    for (int i = 0; i < k.size(); i++) {
//        cout << k[i] << endl;
////        if (i > 0) {
//        check = check ^ k[i];
////        }
//    }

//    cout << scaling * check << endl;

}

