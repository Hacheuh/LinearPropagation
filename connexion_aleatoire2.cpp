//Este programa calcula el MSD de un sistema de N partículas. modified version by HG : version light qui marchent pour des grilles; version 4 : on connecte aleatoirement 4 voisins à chaque pas de temps.
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <vector>
//#include <dos.h>
//Comandos para utilizar los números aleatorios.
extern void seedMT2();
extern void seedMT(unsigned long int);
extern double ranMT(void);
using namespace std;

int N,largeur,Vois=4, freq=2;
double h, t, update;
const int nombreEvenement=4;
int evenementControl=0;
//vector<vector<int> > Adj;
vector<int> defSousGroupes(int i, double ti);
vector<vector<int> > StockSousGroup;
vector<int> Etats;
vector<double> gamma_indiv(double alphaAB, double betaAB, double gammaAB, double omegaAB, double alphaBC, double betaBC, double gammaBC, double omegaBC, double alphaCA, double betaCA, double gammaCA, double omegaCA,  double N, vector<int> etat);
vector<int> Liste_etat();
int ComptageEtat(int etat, vector<int> liste, int quest=0);
double ProbAB(double alpha, double beta, double gamma, double omega, double N, int i);
double ProbBC(double alpha, double beta, double gamma, double omega, double N, int i);
double ProbCA(double alpha, double beta, double gamma, double omega, double N, int i);


//Funciones externas a main() para calcular los rates.
double ProbAB(double alpha, double beta, double gamma, double omega, double N, int i)
{
    vector<int> sousGroupeDef=defSousGroupes(i,t);
    double rateAB = 0.0, nA = ComptageEtat(0,sousGroupeDef,1), nB = ComptageEtat(1, sousGroupeDef,1);
    /*printf("Transition AB\n");
    printf("%d %lf %lf\n", i, nA, nB);*/
    if (nA > 0.0)
    {
        rateAB = ((alpha * pow(nB,beta)) + omega)/pow(nA,gamma);
    }
    else if (nA <= 0.0)
    {
        rateAB = ((alpha * pow(nB,beta)) + omega);
    }
    return rateAB;

}

double ProbBC(double alpha, double beta, double gamma, double omega, double N, int i)
{
    vector<int> sousGroupeDef=defSousGroupes(i,t);
    double rateBC = 0.0, nB = ComptageEtat(1, sousGroupeDef,1), nC = ComptageEtat(2, sousGroupeDef,1);
    /*printf("Transition BC\n");
    printf("%d %lf %lf\n", i, nB, nC);*/
    if (nB > 0.0)
    {
        rateBC = ((alpha * pow(nC,beta)) + omega)/pow(nB,gamma);
    }
    else if (nB <= 0.0)
    {
        rateBC = ((alpha * pow(nC,beta)) + omega);
    }
    return rateBC;

}

double ProbCA(double alpha, double beta, double gamma, double omega, double N, int i)
{
    vector<int> sousGroupeDef=defSousGroupes(i,t);
    double rateCA = 0.0, nC = ComptageEtat(2, sousGroupeDef,1), nA = ComptageEtat(0, sousGroupeDef,1);
    /*printf("Transition CA\n");
    printf("%d %lf %lf\n", i, nC, nA);*/
    if (nC > 0.0)
    {
        rateCA = ((alpha * pow(nA,beta)) + omega)/pow(nC,gamma);
    }
    else if (nC <= 0.0)
    {
        rateCA = ((alpha * pow(nA,beta)) + omega);
    }
    return rateCA;


}


//Se leen los parámetros en la línea de comando.
int main(int argc, const char * argv[])
{
    //Es para utilizar los números aleatorios.
    unsigned int uin = time(NULL);
    srand(uin);
    seedMT(uin);

    int i,j,k,l;
    //int N;
    double T,particiones,omegaAB,omegaBC,omegaCA,alphaAB,alphaBC,alphaCA,betaAB,betaBC,betaCA,gammaAB,gammaBC,gammaCA;

    //Se asignan las entradas de la línea de comando.
    if (argc > 1)
    {
        char * pEnd;
        N           = atoi	(argv[1]);
        T           = strtod(argv[2],&pEnd);
        particiones = strtod(argv[3],&pEnd);
        alphaAB		= strtod(argv[4],&pEnd);
        alphaBC		= strtod(argv[5],&pEnd);
        alphaCA		= strtod(argv[6],&pEnd);
        betaAB		= strtod(argv[7],&pEnd);
        betaBC		= strtod(argv[8],&pEnd);
        betaCA		= strtod(argv[9],&pEnd);
        gammaAB		= strtod(argv[10],&pEnd);
        gammaBC		= strtod(argv[11],&pEnd);
        gammaCA		= strtod(argv[12],&pEnd);
        omegaAB		= strtod(argv[13],&pEnd);
        omegaBC		= strtod(argv[14],&pEnd);
        omegaCA		= strtod(argv[15],&pEnd);
        largeur		= strtod(argv[16],&pEnd);
    }

    char filename[1000];

    //Aquí se concatenan los valores para dar el nombre al archivo .dat
    if (argc > 1)
    {
        char onestr[500];
        sprintf(onestr,"%d",N);
        strcpy(filename,"N_"
              );
        strcat(filename,onestr);

        strcat(filename,".dat");
    }

    //Se definen los strings para dar nombre a los diferentes archivos.
    char filename_nm[1000];
    strcpy(filename_nm,"data_nM_treses");
    strcat(filename_nm,filename);
    ofstream toNM(filename_nm);
    toNM.precision(5);
    toNM.setf(ios::scientific,ios::floatfield);

    char filename_Ind[1000];
    strcpy(filename_Ind,"data_Ind_treses");
    strcat(filename_Ind,filename);
    ofstream toInd(filename_Ind);
    toInd.precision(5);
    toInd.setf(ios::scientific,ios::floatfield);

    char filename_Eve[1000];
    strcpy(filename_Eve,"data_Eve_treses");
    strcat(filename_Eve,filename);
    ofstream toEve(filename_Eve);
    toEve.precision(5);
    toEve.setf(ios::scientific,ios::floatfield);

    char filename_Efinal[1000];
    strcpy(filename_Efinal,"data_Efinal_treses");
    strcat(filename_Efinal,filename);
    ofstream toEfinal(filename_Efinal);
    toEve.precision(5);
    toEve.setf(ios::scientific,ios::floatfield);

    //Las variables de los pasos de tiempo.
    t = 0.0;
    h = 0.1/pow(Vois,betaCA);
    update=0;
    
    cout << "h = " << h << "; Nombre d'itération : " << T/h << endl;

    //Las gammas sirven para calcular el algoritmo de Gillespie.
    //double gamma_1 = 0.0, gamma_2 = 0.0, gamma_3 = 0.0, gamma_tot = 0.0;
    double gamma_tot = 0.0;

    //u sirve para calcular un número aleatorio para calcular el paso de tiempo en Gillespie.
    double u = 0.0;

    //N_double sirve para poder dividir cuando se calcula el MSD o la posición del centro de masa.
    double N_double;
    N_double = N;

    //Variables para calcular los números de ocupación.

    Etats=Liste_etat();
    /*for(i=0; i<N; i++)
    {
        printf("%d ", Etats[i]);
    }

    printf("\n");*/

    int nA = ComptageEtat(0, Etats), nB = ComptageEtat(1, Etats), nC = ComptageEtat(2, Etats);
    //printf("%d %d %d\n", nA, nB, nC);

    for (i=0; i<N; i++)
    {
        toInd << Etats[i] << " ";
        if ((i+1)%largeur==0)
        {
            toInd << endl;
        }
    }
    int pro=0;

    while (t < T)
    {
        if(int((t/h)/int(T/h)*100)!=pro)  // Affihage de la progression
        {
            pro=int((t/h)/int(T/h)*100);
            cout << pro << " /100" << endl;
        }
	//cout << t << endl;
        //Se calculan los parámetros de Gillespie.
        gamma_tot=0.0;
	
        vector<double> gammas = gamma_indiv(alphaAB,betaAB,gammaAB,omegaAB,alphaBC,betaBC,gammaBC,omegaBC,alphaCA,betaCA,gammaCA,omegaCA,N_double,Etats);
	update+=(freq*h);
        for (i=0; i<N; i++)
        {
            gamma_tot+=gammas[i];
        }

        toNM << t << "\t \t" << nA << "\t \t" << nB << "\t \t" << nC << endl;

        double gam = 0.0;

        for(i=0; i<N; i++) 	//Boucle sur la population
        {
            u = ranMT();

            if(u < gammas[i]*h) 	//test de transition
            {
                toEve  << t << "\t \t" << i << "\t \t";
                if (Etats[i]==0)
                {
                    Etats[i]=1;
                    nA--;
                    nB++;
                    toEve << "1" << endl;
                    //printf(" AB");
                }
                else if (Etats[i]==1)
                {
                    Etats[i]=2;
                    nB--;
                    nC++;
                    toEve << "2" << endl;
                    //printf(" BC");
                }

                else if (Etats[i]==2)
                {
                    Etats[i]=0;
                    nC--;
                    nA++;
                    toEve << "0" << endl;
                    //printf(" CA");
                }
                else
                {
                    printf("Erreur dans le choix de transitions");
                }
            }
            //printf(" de l'individu %d.\n", increment);
            /*printf("Nouveaux états : ");
            for(i=0; i<N; i++)
            {
                printf("%d ", Etats[i]);
            }
            printf("\n");*/
            //cout << tiempo_gillespie << endl;
            /*if (evenementControl%10000==0 && evenementControl!=0){
                for (i=0;i<N;i++){
                	toEfinal << Etats[i] << " ";
            	    }
                cout << "Lancement Correlation" << endl;
                system("g++ correlation.cpp -o Corre -lm");
                system("./Corre");
            }*/
            evenementControl++;

        }
        /*if (evenementControl>=nombreEvenement){
            printf("Dépassement de l'événement contrôle\n");
        for (i=0;i<N;i++){
        	toEfinal << Etats[i] << " ";
        	if ((i+1)%largeur==0){
            		toInd << endl;
        	}
            }
            break;
        }*/

        t+=h;

    }

    ////////-----------------------------------------------

    toNM.close();

    toInd.close();

    toEve.close();

    cout << "Final del programa." << endl;
    return 0;
}

vector<int> defSousGroupes(int j, double ti)
{
    
    vector<int> sousGroupe;
    vector<int>::iterator it;
    vector<vector<int> >::iterator itplus;
    int i;

    if(ti==update) 				//si on est à t = n*freq*h on va update les voisinages
    {
        if(j==0)
        {
            StockSousGroup.clear();   //On clear le stockage pour le premier passage
        }
        int tirage, nbrTirage=Vois, verif;
        
        //cout << "id" << "\t \t" << j << endl;
        while(nbrTirage>0)
        {
            verif=1;
            tirage=floor(ranMT()*N);
            //cout << "tirage" << "\t \t" << tirage << endl;
            if (tirage==j)
            {
                verif=0;
            }
            else
            {
                for(it=sousGroupe.begin(); it!=sousGroupe.end(); ++it)
                {
                    if(*it!=tirage)
                    {
                        verif*=1;
                    }
                    else
                    {
                        verif*=0;
                    }
                    //cout << *it << endl;
                }
                if(verif==1)
                {
                    sousGroupe.push_back(tirage);
                    nbrTirage--;
                }
            }
        }
        StockSousGroup.push_back(sousGroupe);	//On stock le sous groupe à la file
	
    }
    else  					// Lorsqu'on est pas à la bonne fréquence on prend juste la valeur du sous groupe stockée
    {	
	for(it=StockSousGroup[j].begin(); it!=StockSousGroup[j].end();++it)
	{
        	sousGroupe.push_back(*it);
	}
    }
    return sousGroupe;
}

vector<double> gamma_indiv(double alphaAB, double betaAB, double gammaAB, double omegaAB, double alphaBC, double betaBC, double gammaBC, double omegaBC, double alphaCA, double betaCA, double gammaCA, double omegaCA,  double N, vector<int> etat)
{
    int i;
    vector<double> gammas;
    for (i=0; i<N; i++)
    {
        //cout << "id_g" << "\t \t" << i << endl;
        if (etat[i]==0)
        {
            gammas.push_back(ProbAB(alphaAB,betaAB,gammaAB,omegaAB,N,i));//printf("OK2\n");
        }
        else if (etat[i]==1)
        {
            gammas.push_back(ProbBC(alphaBC,betaBC,gammaBC,omegaBC,N,i));//printf("OK2\n");
        }
        else if (etat[i]==2)
        {
            gammas.push_back(ProbCA(alphaCA,betaCA,gammaCA,omegaCA,N,i));//printf("OK2\n");
        }
    }
    return gammas;
}

vector<int> Liste_etat()
{
    int i,x,y;
    vector<int> Liste;
    for (i=0; i<N; i++)
    {
        y=floor(i/largeur);
        x=i-y*largeur;
        /*if(y>98){
        	Liste.push_back(0);
        }
        else if(y<0){
        	Liste.push_back(1);
        }
        else{
               	Liste.push_back(2);
           	}*/
        Liste.push_back(floor(ranMT()*2.99999));
    }
    return Liste;
}

int ComptageEtat(int etat, vector<int> liste, int quest)
{
    int i;
    int nombreEtat=0;
    if (quest==0)
    {
        for (i=0; i<liste.size(); i++)
        {
            //printf("%d %d\n",liste[i],etat);
            if (liste[i]==etat)
            {
                //printf("OK\n");
                nombreEtat+=1;
            }
        }
    }
    else if(quest==1)
    {
        for (i=0; i<liste.size(); i++)
        {
            //printf("%d %d\n",liste[i],etat);
            if (Etats[liste[i]]==etat)
            {
                //printf("OK\n");
                nombreEtat+=1;
            }
        }
    }
    else
    {
        printf("Erreur ComptageEtat");
    }
    return nombreEtat;
}














