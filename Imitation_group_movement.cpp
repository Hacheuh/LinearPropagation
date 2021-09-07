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

//fonctions

void Initialisation();
void change_state();
vector<double> gamma_indiv(double alphaAB, double betaAB, double gammaAB, double omegaAB, double alphaBC, double betaBC, double gammaBC, double omegaBC, double alphaCA, double betaCA, double gammaCA, double omegaCA,  double N, double mtx[10000][4]);
double ProbAB(double alpha, double beta, double gamma, double omega, double N, int i);
double ProbBC(double alpha, double beta, double gamma, double omega, double N, int i);
double ProbCA(double alpha, double beta, double gamma, double omega, double N, int i);
double ComptageEtat(int etat, vector<int> liste, int quest=0, int indiv_focus=0); // Répertorier les états du voisinage d'un individu focus grâce à son etat, sa liste de voisinage (moduler par le cut-off)

int Find_head(); // fonction qui trouve l'individu de tête
void Move(double temps, int individu); // Faire bouger les individus actifs
double Dist(int i, int j); // Calculer la distance entre i et j
vector<int> defSousGroupes(int i); // Definir le voisinage de i
void Change_orientation_spont(double temps, double h);


//variables
int perc=2, repetition, rep_max=1, rate=100;
double vitesse=1.0, orientation[3]; // nombre moyen d'individu perçu
int N,largeur, head;
int borne_sup, borne_inf;
double T,particiones,omegaAB,omegaBC,omegaCA,alphaAB,alphaBC,alphaCA,betaAB,betaBC,betaCA,gammaAB,gammaBC,gammaCA;
double Positions[10000][4], Positions_temp[10000][4];
vector <int> Etats; //vestige de précédent codes gardé pour comptage etat

int main(int argc, const char * argv[])
{
    //Es para utilizar los números aleatorios.
    unsigned int uin = time(NULL);
    srand(uin);
    seedMT(uin);

    int i,j,k,l;


    //Se asignan las entradas de la línea de comando.
    //if (argc > 1)
    {
        char * pEnd;
        N           = 100;//atoi	(argv[1]);
        T           = 100;//strtod(argv[2],&pEnd);
        particiones = 0;//strtod(argv[3],&pEnd);
        alphaAB		= 1;//strtod(argv[4],&pEnd);
        alphaBC		= 1;//strtod(argv[5],&pEnd);
        alphaCA		= 1;//strtod(argv[6],&pEnd);
        betaAB		= 1;//strtod(argv[7],&pEnd);
        betaBC		= 1;//strtod(argv[8],&pEnd);
        betaCA		= 1;//strtod(argv[9],&pEnd);
        gammaAB		= 0;//strtod(argv[10],&pEnd);
        gammaBC		= 0;//strtod(argv[11],&pEnd);
        gammaCA		= 0;//strtod(argv[12],&pEnd);
        omegaAB		= 0;//strtod(argv[13],&pEnd);
        omegaBC		= 0;//strtod(argv[14],&pEnd);
        omegaCA		= 0;//strtod(argv[15],&pEnd);
        largeur		= 0;//strtod(argv[16],&pEnd);
    }

    char filename[1000];

    //Aquí se concatenan los valores para dar el nombre al archivo .dat
    if (argc > 1)
    {
        char onestr[500];
        char twostr[500];
        char threestr[500];
        sprintf(onestr,"%d",N);
        sprintf(twostr,"%lf",vitesse);
        sprintf(threestr,"%d",rep_max);
        strcpy(filename,"_N_"
              );
        strcat(filename,onestr);
        strcat(filename,"_V_");
        strcat(filename,twostr);
        strcat(filename,"_rep_");
        strcat(filename,threestr);
        strcat(filename,".dat");
    }

    //Se definen los strings para dar nombre a los diferentes archivos.
    char filename_nm[1000];
    strcpy(filename_nm,"data_nM_IGM");
    strcat(filename_nm,filename);
    ofstream toNM(filename_nm);
    toNM.precision(5);
    toNM.setf(ios::scientific,ios::floatfield);

    char filename_Movers[1000];
    strcpy(filename_Movers,"data_Movers_IGM");
    strcat(filename_Movers,filename);
    ofstream toMovers(filename_Movers);
    toMovers.precision(5);
    toMovers.setf(ios::scientific,ios::floatfield);

    for(repetition=0; repetition<rep_max; repetition++)
    {
        cout << "rep :" << " " << repetition+1 << " /" << rep_max << endl;
	
        //Las variables de los pasos de tiempo.
        double t = 0.0, h = 0.1;//0.1/pow(Rang,betaCA);

        if(repetition==0) 	// on affiche le h que lors de la première répétition
        {
            cout << "h = " << h << endl;
        }

        //Las gammas sirven para calcular el algoritmo de Gillespie.
        double gamma_tot = 0.0;

        //u sirve para calcular un número aleatorio para calcular el paso de tiempo en Gillespie.
        double u = 0.0;

        //N_double sirve para poder dividir cuando se calcula el MSD o la posición del centro de masa.
        double N_double;
        N_double = N;

        Initialisation();
        int nA = int(ComptageEtat(0, Etats)), nB = int(ComptageEtat(1, Etats)), nC = int(ComptageEtat(2, Etats)); // initialisation du compteur d'état



        while (t < T)
        {

            int premier_passage=1;

            vector<double> gammas = gamma_indiv(alphaAB,betaAB,gammaAB,omegaAB,alphaBC,betaBC,gammaBC,omegaBC,alphaCA,betaCA,gammaCA,omegaCA,N_double,Positions); // on génére les taux de transitions de toutes les particules à cet instant. Les distances sont calculer dans cet appel. 
	    
            borne_inf=0;
            borne_sup=N;
            Change_orientation_spont( t, h);
            //Find_head();
            for(i=borne_inf; i<borne_sup; i++)
            {
                Move(h,i);
                if(i!=head)
                {
                    u = ranMT();			//on tire aléatoirement
                    //if (i<20){cout << i << "\t \t" << u << " " << gammas[i]*h << "\t \t" << gammas[i] << "\t \t" << endl;}
                    if(u < gammas[i]*h*10) 		//si le tirage est OK
                    {
                        if(Positions[i][3]==0)
                        {
                            Positions[i][3]=1;
			    Positions[i][2]=orientation[1];
                            nB++;
                            nA--;
                        }
                        else if(Positions[i][3]==1)
                        {

                            Positions[i][3]=2;
			    Positions[i][2]=orientation[2];
                            nC++;
                            nB--;

                        }
                        else if(Positions[i][3]==2)
                        {
                            Positions[i][3]=0;
			    Positions[i][2]=orientation[0];
                            nC++;
                            nB--;
                        }
                        else
                        {
                            cout << "Anomalie d'état détectée" << endl;
                            return 0;
                        }
                    }
                }
            }
            toNM << t << "\t \t" << nA << "\t \t" << nB << "\t \t" << nC << endl;
            for(i=0; i<N; i++)
            {	 
	        toMovers << t << " " << Positions[i][0] << " " << Positions[i][1] << " " << Positions[i][2] << " " << Positions[i][3] << endl;
            }

            t+=h;
	    toMovers << endl;
        }
	
    }
}



//Funciones externas a main() para calcular los rates.

vector<double> gamma_indiv(double alphaAB, double betaAB, double gammaAB, double omegaAB, double alphaBC, double betaBC, double gammaBC, double omegaBC, double alphaCA, double betaCA, double gammaCA, double omegaCA,  double N, double mtx[10000][4])
{   //cout << "lancement gamma" << endl;
    int i;
    vector<double> gammas;
    gammas.clear();
    
    for (i=0; i<N; i++)
    {
        if (mtx[i][3]==0)
        {
            gammas.push_back(ProbAB(alphaAB,betaAB,gammaAB,omegaAB,N,i));//printf("OK2\n");
        }
        else if (mtx[i][3]==1)
        {
            gammas.push_back(ProbBC(alphaBC,betaBC,gammaBC,omegaBC,N,i));//printf("OK2\n");
        }
        else if (mtx[i][3]==2)
        {
            gammas.push_back(ProbCA(alphaCA,betaCA,gammaCA,omegaCA,N,i));//printf("OK2\n");
        }
    }
    return gammas;
}

double ProbAB(double alpha, double beta, double gamma, double omega, double N, int i)
{   //cout << "lncmnt ProbAB" << endl;
    double rateAB = 0.0, nA = ComptageEtat(0, defSousGroupes(i),1,i), nB = ComptageEtat(1, defSousGroupes(i),1,i);

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
{   //cout << "lncmnt ProBC" << endl;
    double rateBC = 0.0, nB = ComptageEtat(1, defSousGroupes(i),1,i), nC = ComptageEtat(2, defSousGroupes(i),1,i);

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
{   //cout << "lncmnt ProbCA" << endl;
    double rateCA = 0.0, nC = ComptageEtat(2, defSousGroupes(i),1,i), nA = ComptageEtat(0, defSousGroupes(i),1,i);

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

vector<int> defSousGroupes(int j)
{   //cout << "lancement defsousgroupes" << endl;
    vector<int> sousGroupe;
    sousGroupe.clear();
    int i;
    double tirage;
    //periodicité horizontale topologique
    /*for(i=1; i<Rang+1; i++) // Attention ce rang doit être un int
    {
        if(j-i>-1)
        {
            sousGroupe.push_back(j-i);
        }
        if (j+i<10000)
        {
            sousGroupe.push_back(j+i);
        }
    }*/
    //cout << endl;
    for(i=0; i<N; i++) //métrique
    {
        if (i!=j)
        {
            //strong cut-off :
            if( Dist(i,j)<=10*perc)
            {
                sousGroupe.push_back(i); //printf("SS grp de %d prend %d\n",j,i);
            }
        }
    }
    return sousGroupe;
}

double ComptageEtat(int etat, vector<int> liste, int quest, int indiv_focus)
{   //cout << "lancement comptage état" << endl;
    int i;
    double nombreEtat=0;
    double tirage;
    if (quest==0)                           //calcul du nombre d'état de manière default, pour l'initialisation des nM
    {
        for (i=0; i<N; i++)
        {

            if (Positions[i][3]==etat)
            {
                nombreEtat+=1;
            }
        }
    }
    else if(quest==1)                      //calcul du nombre d'état de manière controlé, pour le calcul des proba indiv
    {
        //if(Etats[indiv_focus]==2 && Etats[indiv_focus-1]==0){cout << endl;}
        for (i=0; i<liste.size(); i++)
        {
            //if(Etats[indiv_focus]==2 && Etats[indiv_focus-1]==0){printf("%d %d\n",liste[i],Etats[liste[i]]);}
            if (Positions[liste[i]][3]==etat)
            {

                nombreEtat+=exp(-(Dist(indiv_focus,liste[i]))/(perc));
                //if(Etats[indiv_focus]==2 && Etats[indiv_focus-1]==0){cout << nombreEtat << endl;}
            }
        }
    }
    else
    {
        printf("Erreur ComptageEtat");
    }
    //if (indiv_focus==4){cout << indiv_focus << " " << etat << " " << "Nombre état: " << nombreEtat << endl;}
    return nombreEtat;

}



void Move(double temps, int individu)
{
    Positions[individu][0]+=(vitesse*cos(Positions[individu][2])*temps);
    Positions[individu][1]+=(vitesse*sin(Positions[individu][2])*temps);
    //cout << Positions << endl;
    //cout << Positions[0] << endl;
}

double Dist(int i, int j)
{
    double dist;
    if (j==-1)
    {
        dist=sqrt(pow(Positions[i][1],2)+pow(Positions[i][0],2));
    }
    else
    {
        dist=sqrt(pow(Positions[i][1]-Positions[j][1],2)+pow(Positions[i][0]-Positions[j][0],2));
    }

    return dist;
}

void Initialisation()
{ //cout << "lancement initialisation" << endl;
    //Es para utilizar los números aleatorios.
    unsigned int uin = time(NULL);
    srand(uin);
    seedMT(uin);
    int i;
    for (i=0; i<N; i++)
    {
        Positions[i][0]=ranMT()*50;
        Positions[i][1]=ranMT()*50;
        Positions[i][2]=M_PI;
        Positions[i][3]=2;

    }
    orientation[0]=0;
    orientation[1]=0;
    orientation[2]=M_PI;

}

int Find_head()
{//cout << "lancement find_head" << endl;
    int i, j;
    vector<int> liste, liste_comp;
    vector<int>::iterator it, itd, itt;
    liste.clear();
    for(i=0; i<N; i++)
    {
        liste.push_back(i);  // on remplit la liste d'individus qu'on va tester
    }
    for(itd=liste.begin(); itd!=liste.end(); ++itd) // on boucle sur cette liste
    {
        j=0;
        liste_comp.clear();
        for(itt=liste.begin(); itt!=liste.end(); ++itt) //on remplit la liste de comparaison qui est rechargée avec la liste updatée
        {
            liste_comp.push_back(*itt);
        }
        for(it=liste_comp.begin(); it!=liste_comp.end(); ++it) // on boucle sur la liste de comparaison
        {
            if( ( (cos(Positions[*itd][2]) * (Positions[*it][0]-Positions[*itd][0]) + sin(Positions[*itd][2]) * (Positions[*it][1]-Positions[*itd][1]) )/ (Dist(*it,*itd)*vitesse) ) <0 && *itd!=*it)
                    //condition pour que le cosinus entre l'angle directionnel et le vecteur diff de position soit négatif i.e. l'individu  est derrière moi      
            {
		if(*it<*itd){--itd;} // il est important que le pointeur de la première liste soit décrémenté pour compensé l'effaçage 
                liste.erase(liste.begin()+j);  // on élimine l'individu
                j-=1;
            }
            j+=1;
            if (liste.size()<=1) // on sort de la boucle lorsqu'il ne reste plus qu'un individu dans la liste.
            {
                break;
            }
        }
        if (liste.size()<=1)
        {
            break;
        }
    }
    //cout << liste[0] << endl;
    return liste[0];
}

void Change_orientation_spont(double temps, double h)
{//cout << "lancement chgmt orientation" << endl;
    int i, test=1;
    
    if(int(temps/h)%rate==0)
    {	for(i=0;i<N-1;i++){if(Positions[i][3]-Positions[i+1][3]==0){}else{cout << "No go" << endl;test=0;break;}}
	//cout << test << endl;;
	if(test==1){
        head=Find_head();
        //cout << head << " -> changement"<< endl;
        Positions[head][2]=ranMT()*2*M_PI;

        switch (int(Positions[head][3]))
        {
        case 0:
            Positions[head][3]=1;
            break;
        case 1:
            Positions[head][3]=2;
            break;
        case 2:
            Positions[head][3]=0;
            break;
        default:
            cout << "Anomalie dans le changement d'état du leader" << endl;
            break;
        }
    	orientation[int(Positions[head][3])]=Positions[head][2]; // orientation associé à l'état dans lequel passe le leader est changée
    }}
}


