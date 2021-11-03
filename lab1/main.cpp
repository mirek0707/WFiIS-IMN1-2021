#include<iostream>
#include<cmath>
#include<fstream>
using namespace std;

void Euler(double deltat, fstream & f)
{
    const int lambda=-1, tmin=0, tmax=5;
    double y0=1, y_an, y_num;
    for(double t=tmin; t<=tmax; t+=deltat)
    {
        y_an=exp(lambda*t);
        if (t==tmin)
            y_num=y0;
        else
            y_num=y_num+deltat*lambda*y_num;
        f<<t<<"\t"<<y_an<<"\t"<<y_num<<"\t"<<y_num-y_an<<endl;
    }
}

void RK2(double deltat, fstream & f)
{
    const int lambda=-1, tmin=0, tmax=5;
    double y0=1, y_an, y_num, k1, k2;
    for(double t=tmin; t<=tmax; t+=deltat)
    {
        y_an=exp(lambda*t);
        if (t==tmin)
            y_num=y0;
        else
            y_num=y_num+deltat/2*(k1+k2);
        k1=lambda*y_num;
        k2=lambda*(y_num+deltat*k1);
        f<<t<<"\t"<<y_an<<"\t"<<y_num<<"\t"<<y_num-y_an<<endl;
    }
}

void RK4(double deltat, fstream & f)
{
    const int lambda=-1, tmin=0, tmax=5;
    double y0=1, y_an, y_num, k1, k2, k3, k4;
    for(double t=tmin; t<=tmax; t+=deltat)
    {
        y_an=exp(lambda*t);
        if (t==tmin)
            y_num=y0;
        else
            y_num=y_num+deltat*(k1+2*k2+2*k3+k4)/6;
        k1=lambda*y_num;
        k2=lambda*(y_num+deltat*k1/2);
        k3=lambda*(y_num+deltat*k2/2);
        k4=lambda*(y_num+deltat*k3);
        f<<t<<"\t"<<y_an<<"\t"<<y_num<<"\t"<<y_num-y_an<<endl;
    }
}

double V(double t, double omegaV)
{
    return 10*sin(omegaV*t);
}

void RRZ2(double oV, fstream & f)
{
    double deltat=pow(10,-4), L=0.1, C=0.001;
    int R=100;
    double omega0=1/sqrt(L*C);
    double T0=2*M_PI/omega0;
    double tmin=0, tmax=4*T0, omegaV=oV*omega0;
    int Q0=0, I0=0;
    double k1Q, k1I, k2Q, k2I, k3Q, k3I, k4Q, k4I, Q, I;
    for(double t=tmin; t<=tmax; t+=deltat)
    {
        if (t==tmin)
        {
            Q=Q0; I=I0;
            k1Q=I0;
            k1I=V(t,omegaV)/L - Q0/(L*C) - R*I0/L;
            k2Q=I0 + deltat*k1I/2;
            k2I=V(t+0.5*deltat,omegaV)/L - (Q0 + deltat*k1Q/2)/(L*C) - R*(I0 + deltat*k1I/2)/L;
            k3Q=I0 + deltat*k2I/2;
            k3I=V(t+0.5*deltat,omegaV)/L - (Q0 + deltat*k2Q/2)/(L*C) - R*(I0 + deltat*k2I/2)/L;
            k4Q=I0 + deltat*k3I;
            k4I=V(t+deltat,omegaV)/L - (Q0 + deltat*k3Q)/(L*C) - R*(I0 + deltat*k3I)/L;
        }
        else
        {
            Q=Q+deltat*(k1Q+k2Q+k3Q+k4Q)/6;
            I=I+deltat*(k1I+k2I+k3I+k4I)/6;
            k1Q=I;
            k1I=V(t,omegaV)/L - Q/(L*C) - R*I/L;
            k2Q=I + deltat*k1I/2;
            k2I=V(t+0.5*deltat,omegaV)/L - (Q + deltat*k1Q/2)/(L*C) - R*(I + deltat*k1I/2)/L;
            k3Q=I + deltat*k2I/2;
            k3I=V(t+0.5*deltat,omegaV)/L - (Q + deltat*k2Q/2)/(L*C) - R*(I + deltat*k2I/2)/L;
            k4Q=I + deltat*k3I;
            k4I=V(t+deltat,omegaV)/L - (Q + deltat*k3Q)/(L*C) - R*(I + deltat*k3I)/L;
        }
        f<<t<<"\t"<<Q<<"\t"<<I<<endl;
    }
}

int main()
{
    fstream p1, p2, p3, p4;
    p1.open("euler_0.01.txt", ios::out);
    p2.open("euler_0.1.txt", ios::out);
    p3.open("euler_1.txt", ios::out);
    if (p1.good() && p2.good() && p3.good())
    {
        Euler(0.01, p1);
        Euler(0.1, p2);
        Euler(1, p3);
    }
    p1.close();
    p2.close();
    p3.close();

    p1.open("RK2_0.01.txt", ios::out);
    p2.open("RK2_0.1.txt", ios::out);
    p3.open("RK2_1.txt", ios::out);
    if (p1.good() && p2.good() && p3.good())
    {
        RK2(0.01, p1);
        RK2(0.1, p2);
        RK2(1, p3);
    }
    p1.close();
    p2.close();
    p3.close();

    p1.open("RK4_0.01.txt", ios::out);
    p2.open("RK4_0.1.txt", ios::out);
    p3.open("RK4_1.txt", ios::out);
    if (p1.good() && p2.good() && p3.good())
    {
        RK4(0.01, p1);
        RK4(0.1, p2);
        RK4(1, p3);
    }
    p1.close();
    p2.close();
    p3.close();

    p1.open("RRZ_0.5.txt", ios::out);
    p2.open("RRZ_0.8.txt", ios::out);
    p3.open("RRZ_1.txt", ios::out);
    p4.open("RRZ_1.2.txt", ios::out);
    if (p1.good() && p2.good() && p3.good() && p4.good())
    {
        RRZ2(0.5, p1);
        RRZ2(0.8, p2);
        RRZ2(1, p3);
        RRZ2(1.2, p4);
    }
    p1.close();
    p2.close();
    p3.close();
    p4.close();
    return 0;
}