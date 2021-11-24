#include<iostream>
#include<cmath>
#include<fstream>
using namespace std;

const double delta= 0.2;
const int nx = 128; 
const int ny = 128;
const double xmax = delta*nx;
const double ymax = delta*ny;
const double TOL = pow(10,-8);

double V[nx + 1][ny + 1];
int iter = 0;

void relaksacja(fstream &p1, fstream &p2, int k)
{
    double S_act=0, S_prev=0;

    do
    {
        for (int i=k; i<nx; i+=k)
        {
            for(int j=k; j<ny; j+=k)
            {
                V[i][j] = ( V[i+k][j] + V[i-k][j] + V[i][j+k] + V[i][j-k] )/4;
            }
        }
        S_prev=S_act;
        S_act=0;
        for (int i = 0; i <= nx - k; i += k)
        {
            for (int j = 0; j <= ny - k; j += k)
            {
                S_act += pow(k * delta, 2) * (pow((V[i + k][j] - V[i][j]) / (2 * k * delta) + (V[i + k][j + k] - V[i][j + k]) / (2 * k * delta), 2) + pow((V[i][j + k] - V[i][j]) / (2 * k * delta) + (V[i + k][j + k] - V[i + k][j]) / (2 * k * delta), 2)) / 2;
            }
        }
        p1<<iter<<"\t"<<S_act<<endl;
        iter++;
    } while ( abs((S_act-S_prev)/S_prev)>=TOL );

    for (int i = 0; i <= nx; i += k)
    {
        for (int j = 0; j <= ny; j += k)
        {
            p2<<i*delta<<"\t"<<j*delta<<"\t"<<V[i][j]<<endl;
        }
    }
    if(k!=1)
    {
        for (int i = 0; i < nx; i += k)
        {
            for (int j = 0; j < ny; j += k)
            {
                V[i+k/2][j+k/2] = (V[i][j] + V[i+k][j] + V[i][j+k] + V[i+k][j+k] )/4;
                if(i!=nx-k)
                    V[i+k][j+k/2] = (V[i+k][j] + V[i+k][j+k] )/2;
                if (j != ny - k)
                    V[i + k / 2][j + k] = (V[i][j + k] + V[i + k][j + k]) / 2;
                if (j!=0)
                    V[i + k / 2][j] = (V[i][j] + V[i + k][j]) / 2;
                if (i != 0)
                    V[i][j + k / 2] = (V[i][j] + V[i][j + k]) / 2;
            }
        }
    }
}
int main()
{
    for (int x = 0; x <= nx; x++)
    {
        V[x][ny] = (-1) * sin(2 * M_PI * x*delta / xmax);
        V[x][0] = sin(2 * M_PI * x*delta / xmax);
    }
    for (int y = 0; y <= ny; y++)
    {
        V[0][y] = sin(M_PI * y*delta / ymax);
        V[nx][y] = sin(M_PI * y*delta / ymax);
    }

    fstream V16, V8, V4, V2, V1;
    fstream S16, S8, S4, S2, S1;

    V16.open("V_16.txt", ios::out);
    S16.open("S_16.txt", ios::out);
    if (V16.good() &&S16.good())
    {
        relaksacja(S16, V16, 16);
    }

    V8.open("V_8.txt", ios::out);
    S8.open("S_8.txt", ios::out);
    if (V8.good() && S8.good())
    {
        relaksacja(S8, V8, 8);
    }

    V4.open("V_4.txt", ios::out);
    S4.open("S_4.txt", ios::out);
    if (V4.good() && S4.good())
    {
        relaksacja(S4, V4, 4);
    }

    V2.open("V_2.txt", ios::out);
    S2.open("S_2.txt", ios::out);
    if (V2.good() && S2.good())
    {
        relaksacja(S2, V2, 2);
    }

    V1.open("V_1.txt", ios::out);
    S1.open("S_1.txt", ios::out);
    if (V1.good() && S1.good())
    {
        relaksacja(S1, V1, 1);
    }

    V1.close();
    V2.close();
    V4.close();
    V8.close();
    V16.close();

    S1.close();
    S2.close();
    S4.close();
    S8.close();
    S16.close();

    return 0;
}