#include <iostream>
#include <multivector.h>
using namespace std;

int main(int argc, char *argv[]) {

    multivector<double> M = 2 * e(4) - 4 * (e(1));
    cout << M << endl;

}

