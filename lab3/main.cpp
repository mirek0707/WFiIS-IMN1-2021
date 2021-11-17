#include<iostream>
#include<cmath>
#include<fstream>
using namespace std;

const double x0 = 0.01;
const double v0 = 0;
const double delta_t0 = 1;
const double S = 0.75;
const double p = 2; 
const double tmax = 40; 
const double alfa = 5;

struct xv
{
    double x;
    double v;
};

double f(double v)
{
    return v;
}
double g(double x, double v)
{
    return alfa*(1-pow(x,2))*v-x;
}

void kKrokuCzasowego(fstream &plik, double TOL, xv (*metoda)(double, double, double))
{
    double t=0;
    double delta_t=delta_t0;
    double xn=x0;
    double vn=v0;
    double Ex=0;
    double Ev=0;
    xv xv2, xv1;

    plik << t << "\t" << delta_t << "\t" << xn << "\t" << vn << endl;
    do
    {
        xv2=metoda(xn, vn, delta_t);
        xv2=metoda(xv2.x, xv2.v, delta_t);

        xv1=metoda(xn, vn, 2*delta_t);

        Ex=(xv2.x-xv1.x)/(pow(2,p)-1);
        Ev=(xv2.v-xv1.v)/(pow(2,p)-1);

        if (max(abs(Ex),abs(Ev))<TOL)
        {
            t+=2*delta_t;
            xn=xv2.x;
            vn=xv2.v;
            plik << t << "\t" << delta_t << "\t" << xn << "\t" << vn << endl;
        }
        delta_t = pow(((S * TOL) / max(abs(Ex), abs(Ev))),(1/(p+1)))*delta_t;
    }
    while (t<(tmax));
}

double F(double xn1, double xn, double vn1, double vn, double delta_t)
{
    return xn1 - xn - (delta_t/2)*(f(vn)+f(vn1));
}
double G(double xn1, double xn, double vn1, double vn, double delta_t)
{
    return vn1 - vn - (delta_t / 2) * (g(xn,vn) + g(xn1,vn1));
}

double a11=1;
double a12(double delta_t)
{
    return (-1)*delta_t/2;
}
double a21(double delta_t, double x, double v)
{
    return ((-1)* delta_t / 2)*((-2)*alfa*x*v-1);
}
double a22(double delta_t, double x)
{
    return 1-(delta_t/2)*alfa*(1-pow(x,2));
}

xv trapez(double xn, double vn, double delta_t)
{
    xv xv1;
    xv1.x=xn;
    xv1.v=vn;
    double delta=pow(10,-10);
    double dx, dv;
    do
    {
        dx=((-1)*F(xv1.x,xn,xv1.v,vn,delta_t)*a22(delta_t,xv1.x)-(-1)*G(xv1.x,xn,xv1.v,vn,delta_t)*a12(delta_t)) / (a11*a22(delta_t, xv1.x)-a12(delta_t)*a21(delta_t, xv1.x, xv1.v));
        dv=(a11*(-1)*G(xv1.x,xn,xv1.v,vn,delta_t)-a21(delta_t, xv1.x, xv1.v)*(-1)*F(xv1.x,xn,xv1.v,vn,delta_t)) / (a11*a22(delta_t, xv1.x)-a12(delta_t)*a21(delta_t, xv1.x, xv1.v));
        xv1.x+=dx;
        xv1.v+=dv;
    } 
    while (abs(dx)>=delta || abs(dv)>=delta);
    return xv1;    
}

xv RK2(double xn, double vn, double delta_t)
{
    double k1x=f(vn);
    double k1v=g(xn, vn);

    double k2x = f(vn+delta_t*k1v);
    double k2v = g(xn + delta_t * k1x, vn + delta_t * k1v);

    xv xv1;
    xv1.x=xn+(delta_t/2)*(k1x+k2x);
    xv1.v=vn+(delta_t/2)*(k1v+k2v);
    return xv1;
}

    int main()
{
    double TOL1=pow(10,-2), TOL2=pow(10,-5); 
    fstream p1, p2, p3, p4;
    p1.open("trapez1.txt", ios::out);
    if (p1.good())
    {
        kKrokuCzasowego(p1, TOL1, trapez);
    }
    p2.open("trapez2.txt", ios::out);
    if (p2.good())
    {
        kKrokuCzasowego(p2, TOL2, trapez);
    }
    p3.open("RK2_1.txt", ios::out);
    if (p3.good())
    {
        kKrokuCzasowego(p3, TOL1, RK2);
    }
    p4.open("RK2_2.txt", ios::out);
    if (p4.good())
    {
        kKrokuCzasowego(p4, TOL2, RK2);
    }
    p1.close();
    p2.close();
    p3.close();
    p4.close();
    return 0;
}