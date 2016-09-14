#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

void general_tridiag(<vec> &a, <vec> &b, <vec> &c, <vec> &y)

int main(int argc, char *argv[]){
    cout << "Armadillo version: " << arma_version::as_string() << endl;
    /* int version = atoi(argv[1]) */;
    int n = 5;

    double x0 = 0.0;
    double x1 = 1.0;
    double h = (x1 - x0)/(n+1.0);

    float a [n] = {-1.0};
    float b [n+1] = {2.0};
    float c [n] = {-1.0};
}

void general_tridiag(<vec> &a, <vec> &b, <vec> &c, <vec> &y){
    n = len(y);
    u = 0;
}
