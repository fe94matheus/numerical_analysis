#include<vector>
#include "romberg.h"


using std::vector;

void print_matrix(std::vector<std::vector<double>> &A);
void print_vector(vector<double> &A);

vector<double> difereca_vect(vector<double> u, vector<double> v);


double modulo(double x);

double norm_inf(vector<double> u);

vector<double> gauss_banded(vector<vector<double>> &A, vector<double> b, int n, int p);

vector<vector<double>> tridiagonal(vector<double> &d, int m, int p, double a, double b);

vector<double> u(int m, double a, double b);



