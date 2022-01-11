#include <iostream>
#include <cmath>
#include <fstream>
#include "gsl/gsl_math.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_blas.h"

using namespace std;

const int nx = 40;
const int ny = 40;
const int N = (nx + 1)*(ny + 1);
const int delta = 1;
const int delta_t = 1;
const int TA = 40;
const int TB = 0;
const int TC = 30;
const int TD = 0;
double kB = 0.1;
double kD = 0.6;
const int IT_MAX = 2000;

int l(int i, int j)
{
    return i + j*(nx + 1);
} 

int main()
{
    gsl_matrix *A = gsl_matrix_calloc(N, N);
    gsl_matrix *B = gsl_matrix_calloc(N, N);
    gsl_vector *c = gsl_vector_calloc(N);
    gsl_vector *T = gsl_vector_calloc(N);

    //Niezerowe elementy macierzowe z uwzględnieniem WB
    //Wnętrze obszaru (szary obszar na rysunku)
    for (int i = 1; i < nx; i++)
    {
        for (int j = 1; j < ny; j++)
        {
            double a = delta_t / (2 * pow(delta, 2));

            gsl_matrix_set(A, l(i, j), l(i, j) - nx - 1, a);
            gsl_matrix_set(A, l(i, j), l(i, j) - 1, a);
            gsl_matrix_set(A, l(i, j), l(i, j) + 1, a);
            gsl_matrix_set(A, l(i, j), l(i, j) + nx + 1, a);
            gsl_matrix_set(A, l(i, j), l(i, j), -4 * a - 1);

            gsl_matrix_set(B, l(i, j), l(i, j) - nx - 1, -a);
            gsl_matrix_set(B, l(i, j), l(i, j) - 1, -a);
            gsl_matrix_set(B, l(i, j), l(i, j) + 1, -a);
            gsl_matrix_set(B, l(i, j), l(i, j) + nx + 1, -a);
            gsl_matrix_set(B, l(i, j), l(i, j), 4 * a - 1);
        }
    }

    //WB Dirichleta (lewy i prawy brzeg)
    for (int i=0; i<=nx; i+=nx)
    {
        for (int j=0;j<=ny;j++)
        {
            gsl_matrix_set(A, l(i, j), l(i, j), 1);
            gsl_matrix_set(B, l(i, j), l(i, j), 1);
            gsl_vector_set(c, l(i, j), 0);
        }
    }

    //WB von Neumanna na górnym brzegu dla chwili n + 1
    for (int i = 1; i<nx; i++)
    {
        gsl_matrix_set(A, l(i, ny), l(i, ny) - nx - 1, -1/(kB*delta));
        gsl_matrix_set(A, l(i, ny), l(i, ny), 1 + 1/(kB*delta));
        gsl_vector_set(c, l(i, ny), TB);
        for (int j = 0; j <= ny; j++)
        {
            gsl_matrix_set(B, l(i, ny), j, 0);
        }
    }

    //WB von Neumanna na dolnym brzegu dla chwili n + 1
    for (int i = 1; i < nx; i++)
    {
        gsl_matrix_set(A, l(i, 0), l(i, 0), 1 + 1 / (kD * delta));
        gsl_matrix_set(A, l(i, 0), l(i, 0) + nx + 1, -1 / (kD * delta));
        gsl_vector_set(c, l(i, 0), TD);
        for (int j = 0; j <= ny; j++)
        {
            gsl_matrix_set(B, l(i, 0), j, 0);
        }
    }

    //Warunki Początkowe
    for (int j=0;j<=ny;j++)
    {
        gsl_vector_set(T, l(0, j), TA);
        gsl_vector_set(T, l(nx, j), TC);
        for (int i=1; i<nx;i++)
        {
            gsl_vector_set(T, l(i, j), 0);
        }
    }

    //Algorytm CN
    gsl_permutation *p = gsl_permutation_calloc(N);
    int signum;
    gsl_linalg_LU_decomp(A, p, &signum);

    gsl_vector *d = gsl_vector_calloc(N);

    for (int it = 0; it <= IT_MAX; it++)
    {
        gsl_blas_dgemv(CblasNoTrans, 1, B, T, 0, d);
        gsl_blas_daxpy(1, c, d);
        gsl_linalg_LU_solve(A, p, d, T);

        if( it == 100 || it == 200 || it == 500 || it == 1000 || it == 2000 )
        {
            string name = "mapa_T_" + to_string(it) + ".txt";
            fstream p1;
            p1.open(name, ios::out);
            if(p1.good())
            {
                for (int i = 0; i <= nx; i++)
                {
                    for (int j = 0; j <= ny; j++)
                    {
                        p1 << i << "\t" << j << "\t" << gsl_vector_get(T, l(i, j)) << "\n";
                    }
                    p1<<"\n";
                }
            }
            p1.close();

            string name2 = "mapa_rozklad_T_" + to_string(it) + ".txt";
            fstream p2;
            p2.open(name2, ios::out);
            if (p2.good())
            {
                for (int i = 1; i <= nx - 1; i++)
                {
                    for (int j = 1; j <= ny - 1; j++)
                    {
                        p2 << i << "\t" << j << "\t" << ((gsl_vector_get(T, l(i, j) + 1) - 2 * gsl_vector_get(T, l(i, j)) + gsl_vector_get(T, l(i, j) - 1)) / pow(delta, 2)) + ((gsl_vector_get(T, l(i, j) + nx + 1) - 2 * gsl_vector_get(T, l(i, j)) + gsl_vector_get(T, l(i, j) - nx - 1)) / pow(delta, 2)) << "\n";
                    }
                    p2 << "\n";
                }
            }
            p2.close();
        }
    }

    gsl_matrix_free(A);
    gsl_matrix_free(B);
    gsl_vector_free(c);
    gsl_vector_free(T);
    gsl_vector_free(d);
    gsl_permutation_free(p);

    return 0;
}