#include<iostream>
#include<math.h>

double f(double x)
{
	return(8-(4.5*(x-sin(x))));
}

double df(double x)
{
	return(4.5*(cos(x)-1));
}

double newton(double x_0){
    double k,x_1, x_2, x_3,e_0,e_1,e_2;
    int tol=30, i=0;
    
    
    if(f(x_0)==0.0){
        return x_0;
    }
    else{
        i=3;
        while(std::abs(f(x_3))>=1.0e-5&&i<tol){
        x_1=x_0-(f(x_0)/df(x_0));
        x_2=x_1-(f(x_1)/df(x_1));
        x_3=x_2-(f(x_2)/df(x_2));
        e_0=x_1-x_0;
        e_1=x_2-x_1;
        e_2=x_3-x_2;
        k=(log(std::abs(e_2/e_1)))/(log(std::abs(e_1/e_0)));
        x_0=x_3;
        
        
        std::cout<<"k: "<<k<<std::endl;
        
        i++;
        
    }
    std::cout<<"Número de Iterações: "<<i<<std::endl;
    return x_3;

    }

    
}

int main(){
    double r;

    r=newton(3.0);
    std::cout<<"Valor da raiz: "<<r<<std::endl;

}
