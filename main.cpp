#include <iostream>
#include <armadillo>
#include <math.h> //delete if arma can handle exp
#include <time.h>

using namespace std;
using namespace arma;

int writestring2file (char *arg_filename[],
		      char *arg_outstring[],
		      bool arg_delete_file);
int write3vars2file (char *arg_filename[],
		     double *arg_a,
		     double *arg_b,
		     double *arg_c);
void general_tridiag(vec &a, vec &b, vec &c, vec &y);
void specific_tridiag(vec &arg_y);
void LU_decomp(mat &arg_A);

int main(int argc, char *argv[]){
    cout << "Armadillo version: "
	 << arma_version::as_string()
	 << endl;
    /* int version = atoi(argv[1]) */;
    //find cmd-line args
    int n; bool LU
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
      LU = atoi(argv[2]);
      cout << "2 cmd-line args: n="
	   << n << " LU=" << LU
	   << endl;
    }
    else {
      cout << "ERROR you gave to many cmd-line exercise"
	   << endl;
    }
    
    clock_t t0, t1, t2, t3, t_gen, t_spec, t_LU;
    
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
    vec y = h*h*f;
    
    //generate diagonal vector elements
    vec a = ones<vec>(n); a *= -1.0;
    vec b = ones<vec>(n); b *= 2.0;
    vec c = ones<vec>(n); c *= -1.0;
    mat A = eye<mat>(n); A *= 2.0;
    for (i=0; i<n-1; i++){
      A[i,i+1] = -1.0;
      A[i+1,i] = -1.0;
    } //filling A elementwise

    //generate vectors for u
    vec u_gen = zeros<vec>(n);
    vec u_spec = zeros<vec>(n);
    vec u_LU = zeros<vec>(n);

    //allocate parameters for storing data 
    //writestring2file(data_loc, "", 1); //make sure file is deleted
    //writestring2file(data_loc, title, 0); //make new file with new 'title' as first line"
    //writestring2file(data_loc, "h, f2c, f3c", 0);
    char *time_filename[];
    char *u_filename[];
    if (not LU) {
      time_filename = "dderiv_time_c++_nSOMETHING_GS.dat";
      u_filename = "dderiv_u_c++_nSOMETHING_GS.dat";
    }
    else {
      time_filename = "dderiv_time_c++_nSOMETHING_LU.dat";
      u_filename = "dderiv_u_c++_nSOMETHING_LU.dat";
    }
    
    //start exercises
    t0 = clock();
    
    if (not LU) {
      /*calculate u using the general tridiagonal method*/
      //u_gen = general_tridiag(a, b, c, y)
      t1 = clock();
      t_gen = (t1 - t0)/CLOCKS_PER_SEC();
      /*calculate u using the specific tridiagonal method*/
      //u_spec = specific_tridiag(y)
      t2 = clock();
      t_spec = (t2 - t1)/CLOCKS_PER_SEC();
      //write timing-results to file
      //write u(x) to file
    }
    else {
      /*calculate u using LU-decomposition*/
      u_LU = general_LU_decomp(A_matrix=A, vert=y)
      t3 = clock();
      t_LU = (t3 - t0)/CLOCKS_PER_SEC();
      //write timing results to file
      //write u(x) to file
    }

}

void general_tridiag(vec &arg_a, vec &arg_b, vec &arg_c, vec &arg_y){
}
void specific_tridiag(vec &arg_y){
}
void LU_decomp(mat &arg_A){
}

int writestring2file (char *arg_filename[],
		      char *arg_outstring[],
		      bool arg_delete_file) {
  /* Take one single string and add it to the 
     last line of [filename] (defined locally),
     then add endline and close file.
     This is slow and unefficient, but not more is needed.
  */
  if (delete_file) {
    /*delete preexisting file and exit*/
    if( std::remove( filename ) == 0 ) {
      std::cout << "one file successfully deleted:" << std::endl
		<< "\t" << filename << std::endl;
      return 0;
    } //if removed succesfully
    else {
      std::cout << "could not delete file:" << std::endl
		<< "\t" << filename << std::endl;
      return 0;
    } //else not removed succesfully
  } //if boolean 'delete_file'
  
  std::ofstream outfile;
  outfile.open(filename, std::ios::app); 
  outfile << outstring << std::endl;
  outfile.close();
  std::cout << "wrote string to: " << std::endl
	    << "\t" << filename << std::endl;
  return 0;
} // write2file

int write3vars2file (char *arg_filename[],
		     double *arg_a,
		     double *arg_b,
		     double *arg_c) {
  /* Open a file and append three variables to it,
     using the csv-format.*/
  std::ofstream outfile;
  outfile.open(filename, std::ios::app); 
  outfile << std::scientific << std::setprecision(20)
	  << a << ", "
	  << b << ", "
	  << c << std::endl;
  outfile.close();
  /*std::cout << "writing skalars to file: " << std::endl
    << "\t" << filename << std::endl;*/
  return 0;
} // write3vars2file
