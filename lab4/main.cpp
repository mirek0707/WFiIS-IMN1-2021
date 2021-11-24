#include<iostream>
#include<cmath>
#include<fstream>
using namespace std;

const double epsilon = 1;
const double delta = 0.1;
const int nx = 150;
const int ny = 100;
const double V1 = 10;
const double V2 = 0;
const double xmax = delta*nx;
const double ymax = delta*ny;
const double sigmax = 0.1*xmax;
const double sigmay = 0.1*ymax; 
const double TOL = pow(10,-8);

double ro1(double x, double y)
{
    return exp((-1) * (pow((x - 0.35 * xmax) / (sigmax), 2)) - pow((y - 0.5 * ymax) / (sigmay), 2));
}
double ro2(double x, double y)
{
    return -exp((-1) * (pow((x - 0.65 * xmax) / (sigmax), 2)) - pow((y - 0.5 * ymax) / (sigmay), 2));
}

void r_globalna(fstream &Vp, fstream &Sp, fstream &Bp, double wg)
{
    double Vn[nx+1][ny+1] = {0.0};
    double Ro[nx+1][ny+1] = {0.0};
    double Vs[nx+1][ny+1] = {0.0};

    for(int i=0;i<=nx;i++)
    {
        Vn[i][0]=V1;
        Vs[i][0]=V1;
        for (int j=0;j<=ny;j++)
            {
                Ro[i][j] = ro1(i * delta, j * delta) + ro2(i * delta, j * delta);
            }
    }
    int iter=0;
    double S_prev=0.0, S_act=0.0;
    do
    {
        for(int i=1;i<nx;i++)
        {
            for(int j=1;j<ny;j++)
            {
                Vn[i][j]= (Vs[i+1][j] + Vs[i-1][j] + Vs[i][j+1] + Vs[i][j-1] + pow(delta,2)*Ro[i][j]/epsilon)/4;
            }
        }

        for(int j=1;j<ny;j++)
        {
            Vn[0][j]=Vn[1][j];
            Vn[nx][j]=Vn[nx-1][j];
        }

        for (int i = 0; i <= nx; i++)
        {
            for (int j = 0; j <= ny; j++)
            {
                Vs[i][j]=(1-wg)*Vs[i][j]+wg*Vn[i][j];
            }
        }
        S_prev=S_act;
        S_act=0;
        for (int i=0;i<nx;i++)
        {
            for (int j=0;j<ny;j++)
            {
                S_act += pow(delta, 2) * (pow((Vn[i + 1][j] - Vn[i][j]) / delta, 2) / 2 + pow((Vn[i][j+1] - Vn[i][j]) / delta, 2) / 2 - Ro[i][j]*Vn[i][j]);
            }
        }
        Sp<<iter<<"\t"<<S_act<<endl;
        iter++;
    }
    while (abs((S_act - S_prev) / S_prev) >= TOL);

    double Vblad[nx+1][ny+1]={0.0};
    for (int i=1; i<nx; i++)
    {
        for (int j = 1; j < ny; j++)
        {
            Vblad[i][j] = (Vn[i+1][j] - 2*Vn[i][j] + Vn[i-1][j]) / pow(delta,2) + (Vn[i][j+1] - 2*Vn[i][j] + Vn[i][j-1]) / pow(delta,2) + Ro[i][j] / epsilon;
        }
    }
    for (int i=0;i<=nx;i++)
    {
        for (int j=0;j<=ny;j++)
        {
            Vp << i*delta << "\t" << j*delta <<"\t"<<Vn[i][j]<<endl;
            Bp << i*delta << "\t" << j*delta <<"\t"<<Vblad[i][j]<<endl;
        }
    }

}

void r_lokalna(fstream &Sp, double wL)
{
    double V[nx + 1][ny + 1] = {0.0};
    double Ro[nx + 1][ny + 1] = {0.0};

    for (int i = 0; i <= nx; i++)
    {
        V[i][0] = V1;
        for (int j = 0; j <= ny; j++)
        {
            Ro[i][j] = ro1(i * delta, j * delta) + ro2(i * delta, j * delta);
        }
    }
    int iter = 0;
    double S_prev = 0.0, S_act = 0.0;
    do
    {
        for (int i=1;i<nx;i++)
        {
            for(int j=1;j<ny;j++)
            {
                V[i][j]=(1-wL)*V[i][j] + wL*((V[i+1][j] + V[i-1][j] + V[i][j+1] + V[i][j-1] + pow(delta,2)*Ro[i][j]/epsilon))/4;
            }
        }
        for (int j=1; j<ny;j++)
        {
            V[0][j] = V[1][j];
            V[nx][j] = V[nx - 1][j];
        }

        S_prev = S_act;
        S_act = 0.0;
        for (int i = 0; i < nx; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                S_act += pow(delta, 2) * (pow((V[i + 1][j] - V[i][j]) / delta, 2) / 2 + pow((V[i][j + 1] - V[i][j]) / delta, 2) / 2 - Ro[i][j] * V[i][j]);
            }
        }
        Sp << iter << "\t" << S_act << endl;
        iter++;

    } while (abs((S_act - S_prev) / S_prev) >= TOL);
}

int main()
{
    double wg1=0.6, wg2=1.0;
    double wL1 = 1.0, wL2=1.4, wL3=1.8, wL4=1.9;
    fstream p1, p2, p3, p4, p5, p6, p7, p8, p9, p10;
    p1.open("globalna_V_06.txt", ios::out);
    p2.open("globalna_S_06.txt", ios::out);
    p3.open("globalna_blad_06.txt", ios::out);
    if (p1.good()&&p2.good()&&p3.good())
    {
        r_globalna(p1,p2,p3,wg1);
    }

    p4.open("globalna_V_1.txt", ios::out);
    p5.open("globalna_S_1.txt", ios::out);
    p6.open("globalna_blad_1.txt", ios::out);
    if (p4.good() && p5.good() && p6.good())
    {
        r_globalna(p4, p5, p6, wg2);
    }

    p7.open("lokalna_10.txt", ios::out);
    p8.open("lokalna_14.txt", ios::out);
    p9.open("lokalna_18.txt", ios::out);
    p10.open("lokalna_19.txt", ios::out);
    if (p7.good())
    {
        r_lokalna(p7, wL1);
    }
    if (p8.good())
    {
        r_lokalna(p8, wL2);
    }
    if (p9.good())
    {
        r_lokalna(p9, wL3);
    }
    if (p10.good())
    {
        r_lokalna(p10, wL4);
    }

    p1.close();
    p2.close();
    p3.close();
    p4.close();
    p5.close();
    p6.close();
    p7.close();
    p8.close();
    p9.close();
    p10.close();

    return 0;
}