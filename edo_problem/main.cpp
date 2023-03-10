#include<iostream>
#include<cmath>
#include "aproximacao_splines.h"
#include "romberg.h"

//função da base de splines lineares
double phi(double x, double x_0, double x_1, double x_2, double h){
    double y;
    
    if (x<=x_0)
    {
        y=0;
    }
    else if (x>x_0&&x<=x_1)
    {
        y=(x-x_0)/h;
    }
    else if (x>x_1&&x<=x_2)
    {
        y=(x_2-x)/h;
    }
    else if (x>x_2)
    {
        y=0;
    }
    return y;
}

//solução analitica

double u_sol(double x){
    return sin(M_PI*x);
}

/*solução exata para o Teste 2
double u_sol(double x){
    double y;
    double d=0.1;
    if (x>(0.5-d)&&x<(0.5+d))
    {
        y=-(5/203)*x*x+(5/203)*x;
    }
    else{
        y=0.0;
    }
    return y;
}
*/

vector<double> u(int m, double a, double b){
    vector<double> U;
    U.assign(m,0);
    double h=(b-a)/(m+1);

    for (int i = 1; i <=m; i++)
    {
        U[i-1]=u_sol(i*h);
    }
    return U;
}


int main(int argc, char *argv[]){
    vector<vector<double>> A;
    vector<double>c, d , x, U, solc;
    int m, p=1;
    double a, b;
    
    
    m=atoi(argv[1]);
    a=atoi(argv[2]);
    b=atoi(argv[3]);
    double h=(b-a)/(m+1);
    solc.assign(m,0);

    A=tridiagonal(d, m, p, a, b);
    std::cout<<"\nMatriz A tridiagonal\n";
    print_matrix(A);
    
    std::cout<<"\nVetor d\n";
    print_vector(d);
    
    x=gauss_banded(A,d,m,p);
   
    std::cout<<"\nSolução numerica do sistema Ac=d\n";
    print_vector(x);

    std::cout<<"\nSolução exata\n";
    U=u(m,a,b);
    print_vector(U);

    std::cout<<"\nAproximação numerica nos pontos x_i\n";
     
    double solucao;
    for (int j = 1; j <=m; j++)
    {       solucao=0;
            for (int i = 1; i <=m; i++)
            {
                solucao=solucao+x[i-1]*phi(j*h,(i-1)*h,i*h,(i+1)*h,h);
            }
            solc[j-1]=solucao;    
    }
    print_vector(solc);
    std::cout<<"\n";
    
    double erro=norm_inf(difereca_vect(U,solc));
    std::cout<<"\nErro da solução para m="<<m<<": "<<erro<<std::endl;

    return 0;
}