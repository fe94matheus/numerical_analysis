#include<iostream>
#include<vector>
#include<math.h>
#include "romberg.h"

using std::vector;


double romb(double (*f)(double), double a, double b){
    vector<double> T;
    int n, q;
    double soma, h;
    double epsilon=10e-7;
    int ITMAX=30;

    T.assign(ITMAX, 0);
    
    h=(b-a);
    n=1;
    
    T[0]=0.5*h*(f(a)+f(b));
    int k=1;
    do{
        
        soma=0;
        h=0.5*h;
        n=2*n;
        q=1;
        for(int i=1;i<=n-1;i+=2){
            soma+=f(a+i*h);
        }
        T[k]=0.5*T[k-1]+(h*soma);
        for (int i = k-1; i >=0; i--)
        {
            q=q*4;
            T[i]=T[i+1]+((T[i+1]-T[i])/(q-1));
        }
        k++;
        }while(std::abs(T[k-1]-T[k-2])>=epsilon*std::abs(T[k-1])&&k<=ITMAX);
        
        return T[k-1];


}

double romb_2(double (*f)(double, double), double a, double b, double alpha){
    vector<double> T;
    int n, q;
    double soma, h;
    double epsilon=10e-7;
    int ITMAX=30;

    T.assign(ITMAX, 0);
    
    h=(b-a);
    n=1;
    
    T[0]=0.5*h*(f(a,alpha)+f(b,alpha));
    int k=1;
    do{
        
        soma=0;
        h=0.5*h;
        n=2*n;
        q=1;
        for(int i=1;i<=n-1;i+=2){
            soma+=f(a+i*h,alpha);
        }
        T[k]=0.5*T[k-1]+(h*soma);
        for (int i = k-1; i >=0; i--)
        {
            q=q*4;
            T[i]=T[i+1]+((T[i+1]-T[i])/(q-1));
        }
        k++;
        }while(std::abs(T[k-1]-T[k-2])>=epsilon*std::abs(T[k-1])&&k<=ITMAX);
        
        return T[k-1];


}

double romb_3(double (*f)(double, double, double), double a, double b, double alpha, double beta){
    vector<double> T;
    int n, q;
    double soma, h;
    double epsilon=10e-7;
    int ITMAX=30;

    T.assign(ITMAX, 0);
    
    h=(b-a);
    n=1;
    
    T[0]=0.5*h*(f(a,alpha,beta)+f(b,alpha,beta));
    int k=1;
    do{
        
        soma=0;
        h=0.5*h;
        n=2*n;
        q=1;
        for(int i=1;i<=n-1;i+=2){
            soma+=f(a+i*h,alpha,beta);
        }
        T[k]=0.5*T[k-1]+(h*soma);
        for (int i = k-1; i >=0; i--)
        {
            q=q*4;
            T[i]=T[i+1]+((T[i+1]-T[i])/(q-1));
        }
        k++;
        }while(std::abs(T[k-1]-T[k-2])>=epsilon*std::abs(T[k-1])&&k<=ITMAX);
        
        return T[k-1];


}
