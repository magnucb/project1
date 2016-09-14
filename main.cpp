#include <iostream>
#include <armadillo>
#include <math.h>

using namespace std;
using namespace arma;

void general_tridiag(vec &a, vec &b, vec &c, vec &y);
void specific_tridiag(vec &arg_y);
void LU_decomp(mat &arg_A);

int main(int argc, char *argv[]){
    cout << "Armadillo version: " << arma_version::as_string() << endl;
    /* int version = atoi(argv[1]) */;
    int n = 5;

    //generate x-vector
    double x0 = 0.0;
    double x1 = 1.0;
    double h = (x1 - x0)/(n+1.0);
    vec x = linspace<vec>(x0,x1,n);
    vec f = zeros<vec>(n);
    vec f1 = zeros<vec>(n);
    f = 100.0*exp(-10.0*x); //using armadillo
    for (i=0; i<n; i++){
      f1[i] = 100.0*exp(-10.0*x[i]);
    }  //calculating f elementwise
    vec y = h*h*f
    
    //generate diagonal vector elements
    vec a = ones<vec>(n-1); a *= -1.0;
    vec b = ones<vec>(n); b *= 2.0;
    vec c = ones<vec>(n-1); c *= -1.0;
    mat A = eye<mat>(n); A *= 2.0;
    for (int i=0; i<n-1; i++){
      A[i,i+1] = -1.0;
      A[i+1,i] = -1.0;
    } //filling A elementwise
   
}

void general_tridiag(vec &arg_a, vec &arg_b, vec &arg_c, vec &arg_y, int n){
    vec u = zeros<vec>(n);

    for (int i=0; i<n-1; i++){
        k = arg_a[i]/arg_b[i];
        arg_b[i+1] -= k*arg_c[i];
        arg_y[i+1] -= k*arg_y[i];
    }
    for (int i=n-1; i>=0; i--){
        u[i] = (arg_y[i] - u[i+1]*arg_c[i])/arg_b[i];
    }
    return u
}
void specific_tridiag(vec &arg_y, int n){
    vec u = zeros<vec>(n);
    vec d = zeros<vec>(n-2);

    for (int i=1; i<=n-1; i++){
        d[i-1] = (i+1)/i;
    }
    for (int i=1; i<=n-1; i++){
        arg_y[i] -= arg_y[i-1]/d[i-1]
    }
    for (int i=n-1; i>=1; i--){
        u[i] = (arg_y[i] + u[i+1])/d[i]
    }
    return u
}
void LU_decomp(mat &arg_A){
}
/*
x0=0.0; x1=1.0; h=(x1-x0)/(n+1.0)
#vectors of tridiagonal matrix A
a = -1.0*np.ones(n)
b = 2.0*np.ones(n)
c = -1.0*np.ones(n)
A = np.diag(a[1:], -1) + np.diag(b, 0) + np.diag(c[:-1], 1)
#vectors y and x
x = np.linspace(x0,x1,n) 
y = h*h*100*np.exp(-10*x)
u_exact = 1-(1-np.exp(-10))*x-np.exp(-10*x)

t0 = time.clock()

#general gaussian elimination of tri-diagonal matrix
u_gen = general_tridiag(a, b, c, y)
t1 = time.clock()
t_gen = t1 - t0

#specific gaussian elimination of tri-diagonal matrix
u_spec = specific_tridiag(y)
t2 = time.clock()
t_spec = t2 - t1

#LU-decomposition of tri-diagonal matrix
u_LU = general_LU_decomp(A_matrix=A, vert=y)
t3 = time.clock()
t_LU = t3 - t2
*/
