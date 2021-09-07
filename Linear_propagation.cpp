// On regarde avec ce programme la propagation d'un état le long d'une ligne d'agents.


#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <vector>
extern void seedMT2();
extern void seedMT(unsigned long int);
extern double ranMT(void);
using namespace std;

int N,largeur;

const int nombreEvenement=10000;
int evenementControl=1, tester=0;
vector<int> defSousGroupes(int i);
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
    double rateAB = 0.0, nA = ComptageEtat(0, defSousGroupes(i),1), nB = ComptageEtat(1, defSousGroupes(i),1);
    /*printf("Transition AB\n");
    printf("%d %lf %lf\n", i, nA, nB);*/
    if (nA > 0.0){
        rateAB = ((alpha * pow(nB,beta)) + omega)/pow(nA,gamma);
    }
    else if (nA <= 0.0){
        rateAB = ((alpha * pow(nB,beta)) + omega);
    }
    return rateAB;

}

double ProbBC(double alpha, double beta, double gamma, double omega, double N, int i)
{
    double rateBC = 0.0, nB = ComptageEtat(1, defSousGroupes(i),1), nC = ComptageEtat(2, defSousGroupes(i),1);
    /*printf("Transition BC\n");
    printf("%d %lf %lf\n", i, nB, nC);*/
    if (nB > 0.0){
        rateBC = ((alpha * pow(nC,beta)) + omega)/pow(nB,gamma);
    }
    else if (nB <= 0.0){
        rateBC = ((alpha * pow(nC,beta)) + omega);
    }
    return rateBC;

}

double ProbCA(double alpha, double beta, double gamma, double omega, double N, int i)
{
    double rateCA = 0.0, nC = ComptageEtat(2, defSousGroupes(i),1), nA = ComptageEtat(0, defSousGroupes(i),1);
    /*printf("Transition CA\n");
    printf("%d %lf %lf\n", i, nC, nA);*/
    if (nC > 0.0){
        rateCA = ((alpha * pow(nA,beta)) + omega)/pow(nC,gamma);
    }
    else if (nC <= 0.0){
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
    strcpy(filename_nm,"data_nM_lipro");
    strcat(filename_nm,filename);
    ofstream toNM(filename_nm);
    toNM.precision(5);
    toNM.setf(ios::scientific,ios::floatfield);

    char filename_Ind[1000];
    strcpy(filename_Ind,"data_Ind_lipro");
    strcat(filename_Ind,filename);
    ofstream toInd(filename_Ind);
    toInd.precision(5);
    toInd.setf(ios::scientific,ios::floatfield);

    char filename_Eve[1000];
    strcpy(filename_Eve,"data_Eve_lipro");
    strcat(filename_Eve,filename);
    ofstream toEve(filename_Eve);
    toEve.precision(5);
    toEve.setf(ios::scientific,ios::floatfield);
    
    char filename_Efinal[1000];
    strcpy(filename_Efinal,"data_Efinal_lipro");
    strcat(filename_Efinal,filename);
    ofstream toEfinal(filename_Efinal);
    toEve.precision(5);
    toEve.setf(ios::scientific,ios::floatfield);

    char filename_FrontVague[1000];
    strcpy(filename_FrontVague,"data_FrontVague_lipro");
    strcat(filename_FrontVague,filename);
    ofstream toFrontVague(filename_FrontVague);
    toFrontVague.precision(5);
    toFrontVague.setf(ios::scientific,ios::floatfield);


    //Las variables de los pasos de tiempo.
    double 	tiempo_gillespie	= 0.0, dt_gillespie = 0.0;							//Paso de tiempo irregular.
    double 	tiempo_save			= 0.0, dt_save 		= T/particiones;

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
    for (i=0;i<N;i++){
        toInd << Etats[i] << endl;
    }

    while (tiempo_gillespie < T)
    {
	
        //Se calculan los parámetros de Gillespie.
        gamma_tot=0.0;
        vector<double> gammas = gamma_indiv(alphaAB,betaAB,gammaAB,omegaAB,alphaBC,betaBC,gammaBC,omegaBC,alphaCA,betaCA,gammaCA,omegaCA,N_double,Etats);

        for (i=0; i<N; i++)
        {
            gamma_tot+=gammas[i];
        }

        dt_gillespie 	= -1.0*(log(ranMT())/gamma_tot);

        u               = ranMT()*gamma_tot;

        tiempo_gillespie = tiempo_gillespie  + dt_gillespie;

        toNM << tiempo_gillespie << "\t \t" << nA << "\t \t" << nB << "\t \t" << nC << endl;

        toEve << tiempo_gillespie << "\t \t";
	toFrontVague << tiempo_gillespie << "\t \t";

        double gam = 0.0;
        int increment=0;
        while (u > gam)
        {
            gam+=gammas[increment];
            increment+=1;
        }
        toEve << increment-1 << "\t \t";
	if (increment-1>tester){
		tester=increment-1;
	}
	toFrontVague << tester << endl;
        //printf("On choisit la transition");
        if (Etats[increment-1]==0)
        {
            Etats[increment-1]=1;
            nA--;
            nB++;
            toEve << "1" << endl;
            //printf(" AB");
        }
        else if (Etats[increment-1]==1)
        {
            Etats[increment-1]=2;
            nB--;
            nC++;
            toEve << "2" << endl;
            //printf(" BC");
        }

        else if (Etats[increment-1]==2)
        {
            Etats[increment-1]=0;
            nC--;
            nA++;
            toEve << "0" << endl;
            //printf(" CA");
        }
        else{printf("Erreur dans le choix de transitions");}

        evenementControl++;
        if (evenementControl>=nombreEvenement){
            printf("Dépassement de l'événement contrôle\n");
	    for (i=0;i<N;i++){
        	toEfinal << Etats[i] << " ";
        	if ((i+1)%largeur==0){
            		toInd << endl;
        	}
    	    }
            break;
        }

    }

    ////////-----------------------------------------------

    toNM.close();

    toInd.close();

    toEve.close();

    cout << "Fin du programme." << endl;
    return 0;
}

vector<int> defSousGroupes(int j)
{
    vector<int> sousGroupe;
    sousGroupe.clear();
    int i, Rang=1000;
    //periodicité horizontale
    for(i=1;i<Rang+1;i++){
	if (j+i<10000){
    		sousGroupe.push_back(j+i);
	}
	if(j-i>-1){
    		sousGroupe.push_back(j-i);
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
    int i;
    vector<int> Liste;
    for (i=0; i<N; i++)
    {	
	if(i<1){
		Liste.push_back(0);
	}
	else{
	       	Liste.push_back(2);
       	}
    }
    return Liste;
}

int ComptageEtat(int etat, vector<int> liste, int quest)
{
    int i;
    int nombreEtat=0;
    if (quest==0){
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
    else if(quest==1){
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
    else{printf("Erreur ComptageEtat");}
    return nombreEtat;
}














