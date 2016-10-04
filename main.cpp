#include <iostream>
#include <multivector.h>
using namespace std;

int main(int argc, char *argv[]) {

    OrthonormalMetric<double> m;
    multivector<double> M = RCONT(e(1)^e(2), e(1), m);
    cout << M << endl;

}

