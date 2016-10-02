#include <iostream>
#include <multivector.h>
using namespace std;

int main(int argc, char *argv[]) {

    OrthonormalMetric<double> m;
    multivector<double> M = GP(e(1), e(2), m);
    cout << M << endl;

}

