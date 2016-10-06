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
    multivector<double> M = (e(1) ^ (e(2) + 2*e(3)) ^ e(4)) ^ (e(2) + e(3));
    cout << M << endl;

}

