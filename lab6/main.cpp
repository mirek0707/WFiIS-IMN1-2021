#include<iostream>
#include<cmath>
#include<fstream>
#include<algorithm>
#include"mgmres.h"
using namespace std;

const double delta = 0.1;
const int epsilon1 =1;
int epsilon2 = 1;
int V1 = 10;
int V3 = 10;
int V2 =-10;
int V4 =-10;

int N(int nxy)
{
    return (nxy + 1)*(nxy + 1);
}
int j(int nxy, int l)
{
    return floor(l/(nxy+1));
}
int i(int nxy, int l)
{
    return l - j(nxy, l) * (nxy+1);
}

double xymax(int nxy)
{
    return delta*nxy;
}
double sigma(int nxy)
{
    return xymax(nxy)/10.0;
}
double ro1(double x, double y, int nxy)
{
    return exp(-pow((x - 0.25 * xymax(nxy)) / sigma(nxy), 2) - pow((y - 0.5 * xymax(nxy)) / sigma(nxy), 2) );
}
double ro2(double x, double y, int nxy)
{
    return -exp(-pow((x - 0.75 * xymax(nxy)) / sigma(nxy), 2) - pow((y - 0.5 * xymax(nxy)) / sigma(nxy), 2));
} 

double elem( int nxy, int l)
{
    if (i(nxy, l) <= nxy/2)
        return epsilon1;
    else
        return epsilon2;
}

void wypelnianie_D1(double *a, int *ja, int *ia, double *b, int nxy, fstream &p1, fstream &p2)
{
    int k=-1;
    for(int l=0; l<N(nxy); l++)
    {
        int brzeg=0;
        double vb=0;
        if (i(nxy, l)==0)
        {
            brzeg=1;
            vb=V1;
        }
        if (j(nxy, l) == nxy)
        {
            brzeg = 1;
            vb = V2;
        }
        if (i(nxy, l) == nxy)
        {
            brzeg =1;
            vb = V3;
        }
        if (j(nxy, l) == 0)
        {
            brzeg = 1;
            vb = V4;
        }
        //b[l] = -(ro1(delta * i(nxy, l), delta * j(nxy, l), nxy) + ro2(delta * i(nxy, l), delta * j(nxy, l), nxy));
        b[l]=0;
        if (brzeg==1)
        {
            b[l]=vb;
        }
        ia[l] = -1;
        //lewa skrajna przekatna
        if (l-nxy-1>=0 && brzeg==0)
        {
            k++;
            if (ia[l]<0)
                ia[l]=k;
            a[k] = elem(nxy, l) / pow(delta, 2);
            p1 << k<<"\t"<< a[k]<<"\n";
            ja[k] = l - nxy - 1;
        }
        // poddiagonala
        if(l-1>=0 && brzeg==0)
        {
            k++;
            if (ia[l] < 0)
                ia[l] = k;
            a[k] = elem(nxy, l) / pow(delta, 2);
            p1 << k << "\t"<< a[k]<<"\n";
            ja[k] = l - 1;
        }
        //diagonala
        k++;
        if (ia[l] < 0)
            ia[l] = k;
        if (brzeg==0)
        {
            a[k] = -(2 * elem(nxy, l) + elem(nxy, l+1) +elem(nxy, l+nxy+1) ) / pow(delta, 2);
            p1 << k << "\t" << a[k]<<"\n";
        }
        else
        {
            a[k] =1;
            p1 << k << "\t" << a[k]<<"\n";
        }
        ja[k]=l;
        // naddiagonala
        if (l < N(nxy) && brzeg == 0)
        {
            k++;
            a[k] = elem(nxy, l + 1) / pow(delta,2);
            p1 << k << "\t" << a[k] << "\n";
            ja[k]=l+1;
        }
        // prawa skrajna przekątna
        if (l < N(nxy) - nxy - 1 && brzeg == 0)
        {
            k++;
            a[k] = elem(nxy, l + nxy + 1) / pow(delta, 2);
            p1 << k << "\t" << a[k] << "\n";
            ja[k] = l +nxy + 1;
        }
        if (l % 5 == 0 && l != 0)
            p2<<"\n";
        p2 << l << "\t" << i(nxy, l) << "\t" << j(nxy, l) << "\t" << b[l] << "\n";
    }
    int nz_num=k+1;
    ia[N(nxy)]=nz_num;
    for (int z = nz_num; z < 5 * N(nxy); z++)
        p1<<z<<"\t"<<0<<"\n";
}

void met_Poissona1(int nxy, fstream &p1, fstream &p2)
{
    double a[5*N(nxy)];
    int ja[5*N(nxy)];
    int ia[N(nxy)+1];
    fill_n(a, N(nxy) + 1, -1);
    double b[N(nxy)];
    wypelnianie_D1(a, ja, ia, b, nxy, p1, p2);
}

int wypelnianie_D2(double *a, int *ja, int *ia, double *b, int nxy)
{
    int k = -1;
    for (int l = 0; l < N(nxy); l++)
    {
        int brzeg = 0;
        double vb = 0;
        if (i(nxy, l) == 0)
        {
            brzeg = 1;
            vb = V1;
        }
        if (j(nxy, l) == nxy)
        {
            brzeg = 1;
            vb = V2;
        }
        if (i(nxy, l) == nxy)
        {
            brzeg = 1;
            vb = V3;
        }
        if (j(nxy, l) == 0)
        {
            brzeg = 1;
            vb = V4;
        }
        if (V1==V2==V3==V4)
            b[l] = -(ro1(delta * i(nxy, l), delta * j(nxy, l), nxy) + ro2(delta * i(nxy, l), delta * j(nxy, l), nxy));
        else
            b[l] = 0;
        if (brzeg == 1)
        {
            b[l] = vb;
        }
        ia[l] = -1;
        //lewa skrajna przekatna
        if (l - nxy - 1 >= 0 && brzeg == 0)
        {
            k++;
            if (ia[l] < 0)
                ia[l] = k;
            a[k] = elem(nxy, l) / pow(delta, 2);
            ja[k] = l - nxy - 1;
        }
        // poddiagonala
        if (l - 1 >= 0 && brzeg == 0)
        {
            k++;
            if (ia[l] < 0)
                ia[l] = k;
            a[k] = elem(nxy, l) / pow(delta, 2);
            ja[k] = l - 1;
        }
        //diagonala
        k++;
        if (ia[l] < 0)
            ia[l] = k;
        if (brzeg == 0)
        {
            a[k] = -(2 * elem(nxy, l) + elem(nxy, l + 1) + elem(nxy, l + nxy + 1)) / pow(delta, 2);
        }
        else
        {
            a[k] = 1;
        }
        ja[k] = l;
        // naddiagonala
        if (l < N(nxy) && brzeg == 0)
        {
            k++;
            a[k] = elem(nxy, l + 1) / pow(delta, 2);
            ja[k] = l + 1;
        }
        // prawa skrajna przekątna
        if (l < N(nxy) - nxy - 1 && brzeg == 0)
        {
            k++;
            a[k] = elem(nxy, l + nxy + 1) / pow(delta, 2);
            ja[k] = l + nxy + 1;
        }
    }
    int nz_num = k + 1;
    ia[N(nxy)] = nz_num;
    return nz_num;
}

void met_Poissona2(int nxy, fstream &p1)
{
    double a[5 * N(nxy)];
    int ja[5 * N(nxy)];
    int ia[N(nxy) + 1];
    fill_n(a, N(nxy) + 1, -1);
    double b[N(nxy)], V[N(nxy)];
    int nz_num=wypelnianie_D2(a, ja, ia, b, nxy);
    int itr_max = 500;
    int mr = 500;
    double tol_abs = pow(10, -8);
    double tol_rel = pow(10, -8);
    pmgmres_ilu_cr(N(nxy), nz_num, ia, ja, a, V, b, itr_max, mr, tol_abs, tol_rel);
    double max = 0;
    for (int k = 0; k < N(nxy); k++)
    {
        if (delta * i(nxy, k) < max)
        {
            p1<<"\n";
        }
        p1<<delta * i(nxy, k)<<"\t"<<delta * j(nxy, k)<<"\t"<<V[k]<<"\n";
        max = delta * i(nxy, k);
    }
}

int main()
{
    fstream m, w, ma, mb, mc, mp1, mp2, mp3;
    m.open("macierz.txt", ios::out);
    w.open("wektor.txt", ios::out);
    if (m.good() && w.good())
    {
        met_Poissona1(4, m, w);
    }

    ma.open("mapa_a.txt", ios::out);
    if (ma.good())
    {
        met_Poissona2(50, ma);
    }

    mb.open("mapa_b.txt", ios::out);
    if (mb.good())
    {
        met_Poissona2(100, mb);
    }
    mc.open("mapa_c.txt", ios::out);
    if (mc.good())
    {
        met_Poissona2(200, mc);
    }

    V1=V2=V3=V4=0;
    mp1.open("mapa_p1.txt", ios::out);
    if (mp1.good())
    {
        met_Poissona2(100, mp1);
    }

    epsilon2 = 2;
    mp2.open("mapa_p2.txt", ios::out);
    if (mp2.good())
    {
        met_Poissona2(100, mp2);
    }

    epsilon2 = 10;
    mp3.open("mapa_p3.txt", ios::out);
    if (mp3.good())
    {
        met_Poissona2(100, mp3);
    }
    return 0;
}