#ifndef __TOOLS_H
#define __TOOLS_H


#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <math.h>
//#include <gsl/gsl_complex.h>
//#include <gsl/gsl_complex_math.h>
//#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

//#define quickread_DEBUG

using namespace std;

// file array
typedef vector<double> Row;
typedef vector<Row> Table;

//string convertion
std::string toString(double);
std::string toString(int);

// file read and write
bool fexists(string);
Table quickread(string);
int quickread(string,Table);
int quickwrite(string,Table);
// other tools
vector<double> linspace(double,double,int);
vector<double> logspace(double,double,int);
//Table intertable(Table,vector<double>);
//Table intertable(Table,vector<double>,int,int);
void PrintTable(Table);
vector< vector<double> > TableToVector(string);
double FindNearestPointInVector(vector<double>, double);
// additional GSL-like tools
//void gsl_matrix_complex_conjugate(gsl_matrix_complex*);
//void gsl_matrix_complex_print(gsl_matrix_complex*);
//void gsl_matrix_complex_change_basis_UMUC(gsl_matrix_complex*, gsl_matrix_complex*);
//void gsl_matrix_complex_change_basis_UCMU(gsl_matrix_complex*, gsl_matrix_complex*);
// other stuff
double LinInter(double,double, double,double,double);
#endif
