#include<iostream>
#include<vector>
#include<math.h>

using std::vector;


double f(double x){
    //return 1/(1-x);
    //return x*x;
    return 1/(1+(25*x*x));
}



void romb(double a, double b, double epsilon, int ITMAX){
    vector<double> T;
    int n, q;
    double soma, h;

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
        
        std::cout<<"k="<<k-1<<std::endl;

        std::cout<<"T[k]="<<T[k-1]<<std::endl;

}

int main(){

    romb(-1,1,10e-7,30);

    return 0;
}