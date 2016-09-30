#include <iostream>
#include <multivector.h>
using namespace std;

int main(int argc, char *argv[]) {

    multivector<double> M = 1 + (5 * e(2) ^ (3 * e(3)));
    cout << M << endl;

}

