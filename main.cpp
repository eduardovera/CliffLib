#include <iostream>
#include <multivector.h>
using namespace std;

int main(int argc, char *argv[]) {

    multivector<double> M = - 4 * (e(2)) + 3 * e(3) + e(4) + e(5) + e(6) + e(1);
    cout << M << endl;

}

