#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
using namespace std;

const double delta = 0.01;
const int ro = 1;
const int mi = 1;
const int nx = 200;
const int ny = 90;
const int i1 = 50;
const int J1 = 55;
const int j2 = J1+2;
const int IT_MAX = 20000;
const double yJ1 = J1 * delta;
const double yny = ny * delta;

double Qwy(double Qwe)
{
    return Qwe * (pow(yny, 3) - pow(yJ1, 3) - 3 * yJ1 * pow(yny, 2) + 3 * pow(yJ1, 2) * yny) / pow(yny, 3);
}

void WBpsi(vector<vector<double>>& psi, double Qwe)
{
    double y;
    //brzeg A
    for (int j = J1; j <= ny; j++)
    {
        y = j * delta;
        psi[0][j] = Qwe * (pow(y, 3) / 3 - pow(y, 2) * (yJ1 + yny) / 2 + y * yJ1 * yny) / (2 * mi);
    }

    //brzeg C
    for (int j = 0; j <= ny; j++)
    {
        y = j * delta;
        psi[nx][j] = Qwy(Qwe) * (pow(y, 3) / 3 - pow(y, 2) * yny / 2) / (2 * mi) + Qwe * pow(yJ1, 2) * (-yJ1 + 3 * yny) / (12 * mi);
    }
    //brzeg B
    for (int i = 1; i < nx; i++)
    {
        psi[i][ny] = psi[0][ny];
    }
    //brzeg D
    for (int i = i1; i < nx; i++)
    {
        psi[i][0] = psi[0][J1];
    }
    //brzeg E
    for (int j = 1; j <= J1; j++)
    {
        psi[i1][j] = psi[0][J1];
    }
    //brzeg F
    for (int i = 1; i <= i1; i++)
    {
        psi[i][J1] = psi[0][J1];
    }
}
void WBzeta(vector<vector<double>> &zeta, vector<vector<double>> &psi, double Qwe)
{
    double y;
    //brzeg A
    for (int j = J1; j <= ny; j++)
    {
        y = j * delta;
        zeta[0][j] = Qwe * (2 * y - yJ1 - yny) / (2 * mi);
    }
    //brzeg C
    for (int j = 0; j <= ny; j++)
    {
        y = j * delta;
        zeta[nx][j] = Qwy(Qwe) * (2 * y - yny) / (2 * mi);
    }
    //brzeg B
    for (int i = 1; i <= nx - 1; i++)
    {
        zeta[i][ny] = 2 * (psi[i][ny - 1] - psi[i][ny]) / pow(delta, 2);
    }
    //brzeg D
    for (int i = i1 + 1; i <= nx - 1; i++)
    {
        zeta[i][0] = 2 * (psi[i][1] - psi[i][0]) / pow(delta, 2);
    }
    //brzeg E
    for (int j = 1; j <= J1 ; j++)
    {
        zeta[i1][j] = 2 * (psi[i1 + 1][j] - psi[i1][j]) / pow(delta, 2);
    }
    //brzeg F
    for (int i = 1; i <= i1; i++)
    {
        zeta[i][J1] = 2 * (psi[i][J1 + 1] - psi[i][J1]) / pow(delta, 2);
    }
    zeta[i1][J1] = (zeta[i1 - 1][J1] + zeta[i1][J1 - 1]) / 2;
}
void relaksacja(int Qwe, fstream & p1, fstream & p2)
{
    vector<vector<double>> psi(nx+1, vector<double>(ny+1, 0));
    vector<vector<double>> zeta(nx+1, vector<double>(ny+1, 0));

    WBpsi(psi, Qwe);
    int omega;
    for (int it = 0; it < IT_MAX; it++)
    {
        if (it < 2000)
            omega = 0;
        else
            omega = 1;
        for (int i = 1; i < nx ; i++)
        {
            for (int j = 1; j < ny ; j++)
            {
                if (i != 0 && j != 0 && i != nx && j != ny && !(i <= i1 && j == J1) && !(i == i1 && j <= J1) && !(i <= i1 && j < J1))
                {
                    psi[i][j] =  (psi[i + 1][j] + psi[i - 1][j] + psi[i][j + 1] + psi[i][j - 1] - delta * delta * zeta[i][j]) / 4;
                    zeta[i][j] = (zeta[i + 1][j] + zeta[i - 1][j] + zeta[i][j + 1] + zeta[i][j - 1])/4 - omega * ro / (16.0 * mi) * (((psi[i][j + 1] - psi[i][j - 1]) * (zeta[i + 1][j] - zeta[i - 1][j]) - (psi[i + 1][j] - psi[i - 1][j]) * (zeta[i][j + 1] - zeta[i][j - 1])));
                }
            }
        }
        WBzeta(zeta, psi, Qwe);
        double gamma=0;
        for(int i=1; i<nx;i++)
        {
            gamma+=psi[i+1][j2] + psi[i-1][j2] + psi[i][j2+1] + psi[i][j2-1] - 4*psi[i][j2] - pow(delta,2)*zeta[i][j2];
        }
    }
    for (int i = 0; i <= nx; i++)
    {
        for (int j = 0; j <= ny; j++)
        {
            p1<<i * delta<<"\t"<<j * delta<<"\t"<<psi[i][j]<<"\t"<<zeta[i][j]<<"\n";
        }
        p1<<"\n";
    }

    double Vu[nx + 1][ny + 1] = {0};
    double Vv[nx + 1][ny + 1] = {0};

    for (int i = 1; i < nx; i++)
    {
        for (int j = 1; j < ny; j++)
        {
            if (i != 0 && j != 0 && i != nx && j != ny && !(i <= i1 && j == J1) && !(i == i1 && j <= J1) && !(i <= i1 && j < J1))
            {
                Vu[i][j] = (psi[i][j + 1] - psi[i][j - 1]) / (2 * delta);
                Vv[i][j] = -(psi[i+1][j] - psi[i-1][j]) / (2 * delta);
            }
        }
    }
    for (int i = 1; i < nx; i++)
    {
        for (int j = 1; j < ny; j++)
        {
            p2<<i * delta<<"\t"<<j * delta<<"\t"<<Vu[i][j]<<"\t"<<Vv[i][j]<<"\n";
        }
        p2<<"\n";
    }
}

int main()
{
    fstream p1, p2, p3, p4, p5, p6;
    p1.open("relaksacja_-1000.txt", ios::out);
    p2.open("V_-1000.txt", ios::out);
    if(p1.good()&&p2.good())
    {
        relaksacja(-1000, p1, p2);
    }

    p3.open("relaksacja_-4000.txt", ios::out);
    p4.open("V_-4000.txt", ios::out);
    if (p3.good() && p4.good())
    {
        relaksacja(-4000, p3, p4);
    }

    p5.open("relaksacja_4000.txt", ios::out);
    p6.open("V_4000.txt", ios::out);
    if (p5.good() && p6.good())
    {
        relaksacja(4000, p5, p6);
    }
    p1.close();
    p2.close();
    p3.close();
    p4.close();
    p5.close();
    p6.close();
    return 0;
}