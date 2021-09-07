#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <vector>

using namespace std;

const double delta_x = 0.1, x_inf = 0.0, x_sup = 100.0, delta_t = 0.01, t_inf = 0.0, t_sup = 100.0;
const int  x_part = (x_sup-x_inf)/delta_x, t_part = (t_sup-t_inf)/delta_t;
double t;
double eta_I[t_part][x_part], eta_I_old[x_part];
double eta_A[t_part][x_part], eta_A_old[x_part];

double derivative1(int i);
double derivative2(int i);
int perc = 2.0;
double alpha = 1.0, beta = 1.0, v=0.23;


int main ()
{
    cout << "Lancement du programme" << endl;
    ofstream toEtaA("eta_A.dat", ios::out | ios::trunc);
    ofstream toEtaI("eta_I.dat", ios::out | ios::trunc);
    int i, j, t_id;
    if ((delta_t/delta_x)*9.9 < 1)  // condition de convergence
    {

        for(i=int(1.0/delta_x); i<x_part; i++)
        {
            eta_I[0][i]=1.0/x_part; //remplissage uniforme
            eta_I_old[i]=1.0/x_part;
            toEtaI << eta_I[0][i] << "  ";
        }
        toEtaI << endl;
        for(i=0; i<int(1.0/delta_x); i++)
        {
            eta_A[0][i]=1.0/x_part; //remplissage uniforme
            eta_A_old[i]=1.0/x_part;
            toEtaA << eta_A[0][i] << "  ";
        }
        toEtaA << endl;
        for(t=t_inf+delta_t; t<t_sup; t+=delta_t)
        {
            for(i=0; i<x_part; i++)
            {
                t_id=t/delta_t;
                eta_I[t_id][i]=eta_I_old[i]+derivative1(i)*delta_t;
                eta_A[t_id][i]=eta_A_old[i]+derivative2(i)*delta_t;
                toEtaI << eta_I[t_id][i] << "  ";
                toEtaA << eta_A[t_id][i] << "  ";
            }
            for(i=0; i<x_part; i++)
            {
                t_id=t/delta_t;
                eta_I_old[i]=eta_I[t_id][i];
                eta_A_old[i]=eta_A[t_id][i];
            }
            toEtaI << endl;
            toEtaA << endl;
        }
    }
    else
    {
        cout << "Condition de convergence rejetée" << endl;
    }
    cout << "Fin du programme" << endl;
}

double derivative1(int i)
{
    int j, borne_inf, borne_sup;
    double eta_sum = 0.0, f;
    if(i+perc/delta_x>x_part)
    {
        borne_sup=x_part;
    }
    else
    {
        borne_sup=i+perc/delta_x;
    }

    if(i-perc/delta_x<0)
    {
        borne_inf=0;
    }
    else
    {
        borne_inf=i-perc/delta_x;
    }
    for(j=i+1; j<borne_sup+1; j++)
    {
        eta_sum+=(eta_A_old[j]*delta_x);
    }
    for(j=borne_inf; j<i; j++)
    {
        eta_sum+=(eta_A_old[j]*delta_x);
    }
    f=-alpha*pow(eta_sum,beta)*eta_I_old[i];
    return f;
}

double derivative2(int i)
{
    int j, borne_inf, borne_sup;
    double eta_sum = 0.0, f, deriv_etaA_x;
    if(i+perc/delta_x>x_part)
    {
        borne_sup=x_part;
    }
    else
    {
        borne_sup=i+perc/delta_x;
    }

    if(i-perc/delta_x<0)
    {
        borne_inf=0;
    }
    else
    {
        borne_inf=i-perc/delta_x;
    }
    for(j=i+1; j<borne_sup+1; j++)
    {
        eta_sum+=(eta_A_old[j]*delta_x);
    }
    for(j=borne_inf; j<i; j++)
    {
        eta_sum+=(eta_A_old[j]*delta_x);
    }
    switch(i)
    {
    case 0:
        deriv_etaA_x=0;
        break;
    case x_part:
        deriv_etaA_x=0;
        break;
    default:
        deriv_etaA_x=(eta_A_old[i+1]-eta_A_old[i-1])/(2*delta_x);
        break;
    }

    f=alpha*pow(eta_sum,beta)*eta_I_old[i]-v*deriv_etaA_x;
    return f;
}
