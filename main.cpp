#include <iostream>
#include <multivector.h>
using namespace std;

int main(int argc, char *argv[]) {

    multivector<double> M = 5 + e(1) ^ e(2);
    cout << M << endl;

}

