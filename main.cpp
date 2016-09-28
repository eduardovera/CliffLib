#include <iostream>
#include <multivector.h>
using namespace std;

int main(int argc, char *argv[]) {

    multivector<double> M = e(5) ^ e(1) ^ e(2) ^ e(3) ^ e(4) ^ e(6);
    cout << M << endl;


}

