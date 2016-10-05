#include <iostream>
#include <multivector.h>
using namespace std;

int main(int argc, char *argv[]) {

    /**
     *
     * TO DO:
     *
     * - Assert if it's a blade before some operations (e.g. REVERSE, SQR_NORM_REVERSE).
     * - Change pointers to direct access
     * - Make orthogonal metric accept (P,Q,R) pattern
     *
     **/

    OrthonormalMetric<double> m;
    multivector<double> B = e(1);
    multivector<double> C = e(1)^e(2);
    multivector<double> A = GP(B, C, m);
    B = IGP(A, C, m);

    cout << B << endl;

}

