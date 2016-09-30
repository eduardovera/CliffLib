#include <iostream>
#include <multivector.h>
using namespace std;

int main(int argc, char *argv[]) {

    multivector<double> M = (e(2) ^ e(1));
    cout << take_grade(M) << endl;

}

