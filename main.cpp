#include <iostream>
#include <armadillo>
//#include <math.h> //delete if arma can handle exp
#include <time.h>
#include <string>

using namespace std;
using namespace arma;

void writestring2file(string arg_filename, string arg_outstring);
int write3vars2file (string arg_filename, double *arg_a, double *arg_b, double *arg_c);
double general_tridiag(vec &arg_a, vec &arg_b, vec &arg_c, vec &arg_u, vec &arg_y, int n);
double specific_tridiag(vec &arg_u, vec &arg_y, int n);
double LU_decomp(mat &arg_A, vec &arg_y);

int main(int argc, char *argv[]){
    cout << "Armadillo version: "
	 << arma_version::as_string()
	 << endl;
    /* int version = atoi(argv[1]) */;
    //find cmd-line args
    int n; bool LU;
    if (argc == 1){
      n = 5;
      LU = false;
      cout << "No cmd-line args: n="
	   << n << " LU=" << LU
	   << endl;
    }
    else if (argc == 2) {
      n = atoi(argv[1]);
      LU = false;
      cout << "1 cmd-line args: n="
	   << n << " LU=" << LU
	   << endl;
    }
    else if (argc == 3) {
      n = atoi(argv[1]);
      string cmd_arg2(argv[2]);
      if (cmd_arg2 == "1") {
          LU = true;
      };
      cout << "2 cmd-line args: n="
	   << n << " LU=" << LU
	   << endl;
    }
    else {
      cout << "ERROR you gave to many cmd-line exercise"
	   << endl;
    }
    
    double t_gen, t_spec, t_LU; //CPU time in seconds it takes to calculate the diff. methods
    
    //generate x,f,y -vectors
    double x0 = 0.0;
    double x1 = 1.0;
    double h = (x1 - x0)/(n+1.0);
    vec x = linspace<vec>(x0,x1,n);
    vec f = zeros<vec>(n);
    vec f1 = zeros<vec>(n);
    f = 100.0*exp(-10.0*x); //using armadillo
    for (int i=0; i<n; i++){
      f1(i) = 100.0*exp(-10.0*x(i));
    }  //calculating f elementwise
    vec y = h*h*f;

    //generate diagonal vector elements and matrix
    vec a = ones<vec>(n-1); a *= -1.0;
    vec b = ones<vec>(n); b *= 2.0;
    vec c = ones<vec>(n-1); c *= -1.0;
    mat A (n,n, fill::eye); A *= 2.0;
    for (int i=0; i<n-1; i++){
      A(i,i+1) = -1.0;
      A(i+1,i) = -1.0;
    } //filling A elementwise

    //generate vectors for u
    vec u_gen = zeros<vec>(n);
    vec u_spec = zeros<vec>(n);
    vec u_LU = y;

    //allocate parameters for storing data 
<<<<<<< HEAD
    string time_filename = "/data/dderiv_time_c++.dat";
    string u_filename = "/data/dderiv_u_c++_n" + to_string(n) + "_";

    //start exercises
=======
<<<<<<< HEAD
    //writestring2file(data_loc, "", 1); //make sure file is deleted
    //writestring2file(data_loc, title, 0); //make new file with new 'title' as first line"
    //writestring2file(data_loc, "h, f2c, f3c", 0);
    string prename;
    string u_filename; string time_filename; //filenames of two datafiles
    if (LU) {
        time_filename = "dderiv_time_c++_nSOMETHING_LU.dat";
        u_filename = "dderiv_u_c++_nSOMETHING_LU.dat";
    }
    else {
        time_filename = "dderiv_time_c++_n" + "_GS.dat";
        u_filename = "dderiv_u_c++_n" + "_GS.dat";
=======
    string time_filename = "dderiv_time_c++_n" + to_string(n) + "_";
    string u_filename = "dderiv_u_c++_n" + to_string(n) + "_";
>>>>>>> 0a0cf71a84106276c3a5132bbc418dfc7febc889
    if (LU) {
        /*calculate u using LU-decomposition*/
        t_LU = LU_decomp(A, u_LU); // turns empty array U_LU into solution

        //write timing results to file
        string string_time_data = "method: n=" + to_string(n) + " time=";
        writestring2file(time_filename, "LU" + string_time_data + to_string(t_LU));
        //write u(x)-results to file
        u_filename += "LU.dat";
        //remove old files
        //ofstream::outfile.open(u_filename)
        string string_u_data = "x, u_LU";
        writestring2file(u_filename, string_u_data);
        for (int i=0; i < n; i++) {
            string_u_data = to_string(x(i)) + ", " + to_string(u_gen(i)) + ", " + to_string(u_spec(i));
            writestring2file(u_filename, string_u_data);
        }
    } else {
        /*calculate u using the general tridiagonal method*/
        t_gen = general_tridiag(a, b, c, u_gen, y, n); //turns empty array u_gen into solution

        /*calculate u using the specific tridiagonal method*/
        t_spec = specific_tridiag(u_spec, y, n); //turns empty array u_spec into solution

        //TEST
        cout << "TEST:  calculated specific_tridiag" << endl;

        //write timing-results to file
        string string_time_data = "tridiagonal method: n=" + to_string(n) + " time=";
        writestring2file(time_filename, "general" + string_time_data + to_string(t_gen));
        writestring2file(time_filename, "specific" + string_time_data + to_string(t_spec));

        //TEST
        cout << "TEST:  wrote time to datafile" << endl;

        //write u(x)-results to file
        u_filename += "tridiag.dat";
<<<<<<< HEAD
        string string_u_data = "x, u_gen, u_spec";
        ofstream outfile; outfile.open(u_filename.c_str()); outfile.close();

        writestring2file(u_filename, string_u_data);
        for (int i=0; i < n; i++) {
            string_u_data = to_string(x(i)) + ", " + to_string(u_gen(i)) + ", " + to_string(u_spec(i));
            writestring2file(u_filename, string_u_data);
        }
        //TEST
        cout << "TEST:  wrote u to data file" << endl;
=======
>>>>>>> origin/master
    }
    
    //start exercises
    t0 = clock();
    
    if (not LU) {
      /*calculate u using the general tridiagonal method*/
      general_tridiag(a, b, c, u_gen, y, n); //turns empty array u_gen into solution
      t1 = clock();
      t_gen = (t1 - t0)/CLOCKS_PER_SEC;
      /*calculate u using the specific tridiagonal method*/
      specific_tridiag(u_spec, u_spec, y, n); //turns empty array u_spec into solution
      t2 = clock();
      t_spec = (t2 - t1)/CLOCKS_PER_SEC;
      //write timing-results to file
      //write u(x) to file
    }
    else {
      /*calculate u using LU-decomposition*/
      LU_decomp(A, u_LU); // turns empty array U_LU into solution
      t3 = clock();
      t_LU = (t3 - t0)/CLOCKS_PER_SEC;
      //write timing results to file
      //write u(x) to file
>>>>>>> 0a0cf71a84106276c3a5132bbc418dfc7febc889
    }
}

double general_tridiag(vec &arg_a, vec &arg_b, vec &arg_c, vec &arg_u, vec &arg_y, int n){
    double k;
    clock_t t0, t1;

    for (int i=0; i<n-1; i++){
        k = arg_a(i)/( (double) arg_b(i) );
        arg_b(i+1) -= k*arg_c(i);
        arg_y(i+1) -= k*arg_y(i);
    }
    for (int i=n-2; i>=0; i--){
        arg_u(i) = (arg_y(i) - arg_u(i+1)*arg_c(i))/( (double) arg_b(i) );
    }
    return (t1 - t0)/((double) CLOCKS_PER_SEC);
}

double specific_tridiag(vec &arg_u, vec &arg_y, int n){
    //TEST
    cout << "TEST:  started specific tridiag" << endl;

    clock_t t0,t1;
    vec d = zeros<vec>(n);

    //TEST
    cout << "TEST:  started forward subst." << endl;

    for (int i=1; i<=n-1; i++){
        d(i-1) = (i+1)/( (double) i );
        arg_y(i) -= arg_y(i-1)/d(i-1);
        cout << i << endl
             << d(i-1) << endl
             << arg_y(i) << endl;
    }

    //TEST
    cout << "TEST:  ended forward subst." << endl;
    //TEST
    cout << "TEST:  started backward subst." << endl;

    for (int i=n-2; i>=1; i--){
        arg_u(i) = (arg_y(i) + arg_u(i+1))/d(i);
    }

    //TEST
    cout << "TEST:  ended backward subst." << endl;

    return (t1 - t0)/((double) CLOCKS_PER_SEC);
}

double LU_decomp(mat &arg_A, vec &arg_y){
    /*Solve the equation A*u = y were the matrix A
     * is fetched as 'arg_a', and y is fetched as 'arg_y'.
    */
    clock_t t0, t1;

    mat L,U(size(arg_A)); //initialize the matrices required for LU-decomp.
    lu(L,U, arg_A); //calculate lower and upper triangular matrices from A
    solve(L, arg_y); //solve L*w = y for w (where w is stored in the y-array)
    solve(U, arg_y); //solve U*u = w (where u is stored in the w-array(which is stored in the y-array))
    //array of argument 'arg_y' has now become the solution u of 'A*u = y'

    return (t1 - t0)/((double) CLOCKS_PER_SEC);
}

void writestring2file (string arg_filename, string arg_outstring) {
  /* Take one single string and add it to the 
     last line of [filename] (defined locally),
     then add endline and close file.
     This is slow and unefficient, but not more is needed.
  */
  
  ofstream outfile;
  outfile.open(arg_filename, ios::app);
  outfile << arg_outstring << endl;
  outfile.close();
  cout << "wrote string to: " << endl
        << "\t" << arg_filename << endl;
} // writestring2file


/*
void write3vars2file (string arg_filename, double *arg_a, double *arg_b, double *arg_c) {
  /* Open a file and append three variables to it,
     using the csv-format.//
  std::ofstream outfile;
  outfile.open(arg_filename, std::ios::app);
  outfile << std::scientific << std::setprecision(20)
      << arg_a << ", "
      << arg_b << ", "
      << arg_c << std::endl;
  outfile.close();
    } // write3vars2file
*/
