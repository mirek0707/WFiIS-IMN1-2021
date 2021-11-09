#include<iostream>
#include<fstream>
#include<cmath>
using namespace std;
const double beta = 0.001;
const int N = 500;
const double gamma1 = 0.1;
const int tmax = 100;
const double delta_t = 0.1;
const double u0 = 1;
const double TOL = pow(10,-6);
double f(double u)
{
    return (beta*N-gamma1)*u -beta*pow(u,2);
}
double alfa()
{
    return beta*N-gamma1;
}
void Picard(fstream & p)
{
    double un=u0, un1=u0;
    for(double t=0; t<=tmax; t+=delta_t)
    {
        un=un1;
        int i=0;
        while(abs(un1-un)<=TOL && i<20)
        {    
            un1=un + delta_t/2*(f(un)+f(un1));
            i++;
        }
        p<<t<<"\t"<<un<<"\t"<<N-un<<endl;
    }
}
void Newton(fstream & p)
{
    double un=u0, un1=u0;
    for(double t=0; t<=tmax; t+=delta_t)
    {
        un=un1;
        int i=0;
        while(abs(un1-un)<=TOL && i<20)
        {    
            un1=un1 - (un1 - un - delta_t/2*(f(un)+f(un1)))/(1-delta_t/2*(alfa() -2*beta*un1 ));
            i++;
        }
        p<<t<<"\t"<<un<<"\t"<<N-un<<endl;
    }
}
const double a11=1/4;
const double a22=1/4;
const double a12=1/4 - sqrt(3)/6;
const double a21=1/4 + sqrt(3)/6;
double m11(double U1)
{
   return 1-delta_t*a11*(alfa()-2*beta*U1); 
}
double m12(double U2)
{
   return -delta_t*a12*(alfa()-2*beta*U2); 
}
double m21(double U1)
{
   return -delta_t*a21*(alfa()-2*beta*U1); 
}
double m22(double U2)
{
   return 1-delta_t*a22*(alfa()-2*beta*U2); 
}
double F1(double U1, double U2, double un)
{
    return U1 - un + delta_t*(a11*(alfa()*U1 - beta*pow(U1,2)) + a12*(alfa()*U2 - beta*pow(U2,2)));
}
double F2(double U1, double U2, double un)
{
    return U2 - un + delta_t*(a21*(alfa()*U1 - beta*pow(U1,2)) + a22*(alfa()*U2 - beta*pow(U2,2)));
}
void RK2(fstream & p)
{
    double un=u0, un1=u0;
    double U1, U2, delta_U1, delta_U2;
    for(double t=0; t<=tmax; t+=delta_t)
    {
        un=un1;
        int i=0;
        U1=U2=un;
        while(abs(un1-un)<=TOL && i<20)
        {
            delta_U1=(F2(U1,U2,un)*m12(U2) - F1(U1,U2,un)*m22(U2))/(m11(U1)*m22(U2) - m12(U2)*m21(U1));
            delta_U2=(F1(U1,U2,un)*m21(U1) - F2(U1,U2,un)*m11(U1))/(m11(U1)*m22(U2) - m12(U2)*m21(U1));
            U1+=delta_U1;
            U2+=delta_U2;
            un1=un+delta_t/2*( f(U1) +f(U2));
            i++;
        }
        p<<t<<"\t"<<un<<"\t"<<N-un<<endl;
    }

}
int main()
{
    fstream p1, p2, p3;
    p1.open("Picard.txt", ios::out);
    if (p1.good())
    {
        Picard(p1);
    }
    p2.open("Newton.txt", ios::out);
    if (p2.good())
    {
        Newton(p2);
    }
    p3.open("RK2.txt", ios::out);
    if (p3.good())
    {
        RK2(p3);
    }
    p1.close();
    p2.close();
    p3.close();
    return 0;
}