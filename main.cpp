#include <iostream>
#include <multivector.h>
using namespace std;

int main(int argc, char *argv[]) {

    multivector<double> M({0.0, -3.0, 6.0, 5.0, 4.0, 9.0, 5.0, 2.3, 2.1});
    cout << M << endl;

//    cout << e(1) << endl;

//    cout << M.pop(2) << endl;

}

