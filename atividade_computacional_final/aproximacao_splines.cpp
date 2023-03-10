#include<iostream>
#include<vector>
#include<cmath>
#include "aproximacao_splines.h"
#include "romberg.h"

/*funções a serem integradas*/
double g(double x){
    return 1.0;
}


double q(double x){
    return M_PI*M_PI;
}

double f(double x){
    return 2.0*M_PI*M_PI*sin(M_PI*x);
}

/*f para os Testes 3 e 4
double f(double x){
    double f_0;
    double d=5;
    if (x>(0.5-d)&&x<(0.5+d))
    {
        f_0=0.1;
    }
    else{
        f_0=0.0;
    }
    return f_0;
}*/
/*f para o Teste 1
double f(double x){
    return 12.0*x*(1.0-x)-2.0;
}
*/

double func_int_q_sup_inf(double x, double x_alpha, double x_beta){
    return (x_beta-x)*(x-x_alpha)*q(x);
}


double integral_q_pos(double x, double x_alpha){
    return (x-x_alpha)*(x-x_alpha)*q(x);
}

double integral_q_neg(double x, double x_alpha){
    return (x_alpha-x)*(x_alpha-x)*q(x); 
}

double integral_f_pos(double x, double x_alpha){
    return (x-x_alpha)*f(x);
}

double integral_f_neg(double x, double x_alpha){
    return (x_alpha-x)*f(x); 
}



/*-----------------------------------------------------------------------------------*/
vector<double> difereca_vect(vector<double> u, vector<double> v){
    vector<double> w;
    w.assign(u.size(), 0);
    for(int i=0; i<u.size();i++){
        w[i]=u[i]-v[i];
    }
    return w;
}

double modulo(double x){
    if(x<0){
        return -x;
    }
    else{
        return x;
    }
}

double norm_inf(vector<double> u){
    double maior;
    maior=-INT16_MAX;
    for(double x: u){
        if(modulo(x)>maior){
            maior=modulo(x);
        }

    }
    return maior;
}

vector<double> gauss_banded(vector<vector<double>> &A, vector<double> b, int n, int p){
    vector<double> x;
    double m;
    x.assign(n,0);

    for (int k = 0; k < n-1; k++)
    {
        for (int i = 1; i < std::min(p+1,n-k); i++)
        {
            m=A[k+i][p-i]/A[k][p];
            b[k+i]=b[k+i]-m*b[k];
            for (int j = p-i+1; j < 2*p; j++)
            {
                A[k+i][j]=A[k+i][j]-m*A[k][j+1];
            }   
        }   
    }
    
    double sum;
    x[n-1]=b[n-1]/A[n-1][p];
    
    for (int i = n-2; i >=0; i--)
    {
        sum=0;
        for (int j = p+1; j < 2*p+1; j++)
        {
            sum=sum+A[i][j]*x[i+1];
        }
        
        x[i]=(b[i]-sum)/A[i][p];
    }
    
    return x;
}

void print_matrix(vector<vector<double>> &A){
    for (vector<double> linha:A)
    {
        for (double a:linha)
        {
            std::cout<<a<<"\t";
        }
        std::cout<<"\n";
    }
}

void print_vector(vector<double> &A){
    
    for (double a: A)
    {
        std::cout<<a<<"\t";
    }
    std::cout<<"\n";
    
}

vector<vector<double>> tridiagonal(vector<double> &d, int m, int p, double a, double b){
    vector<vector<double>> A(m,vector<double>(2*p+1));
    vector<double> x;
    d.assign(m,0);
    x.assign(m+2,0);
    double h=(b-a)/(m+1);
    
    //pontos x_i's no intervalo [a, b]
    for (int i = 0; i <= m+1; i++)
    {
        x[i]=i*h;
    }
    printf("\n--------pontos x_i----------\n");
    print_vector(x);
    printf("\n------------------------\n");
//A[i][p] elementos da diagonal principal e d[i]'s <u,phi_i> 
    for (int i = 1; i <=m; i++)
    {
        A[i-1][p]=(1/(h*h))*romb(g,x[i-1],x[i])+(1/(h*h))*romb(g,x[i],x[i+1])+
        (1/(h*h))*romb_2(integral_q_pos,x[i-1],x[i],x[i-1])+
        (1/(h*h))*romb_2(integral_q_neg,x[i],x[i+1],x[i+1]);
        d[i-1]=(1/h)*romb_2(integral_f_pos,x[i-1],x[i],x[i-1])+(1/h)*romb_2(integral_f_neg,x[i],x[i+1],x[i+1]);
    }
    //A[i][2*p] elementos da diagonal superior
    for (int i = 1; i <m; i++)
    {
        A[i-1][2*p]=-(1/(h*h))*romb(g,x[i],x[i+1])+(1/(h*h))*romb_3(func_int_q_sup_inf,x[i],x[i+1],x[i],x[i+1]);
    }
    //A[i][0] elementos da diagonal inferior
    for (int i = 2; i <=m; i++)
    {
        A[i-1][0]=-(1/(h*h))*romb(g,x[i-1],x[i])+(1/(h*h))*romb_3(func_int_q_sup_inf,x[i-1],x[i],x[i-1],x[i]);

    }

    return A;
}



