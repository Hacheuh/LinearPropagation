//Première version du script IGM qui montre un groupe en mouvement et une propagation par vague de l'info de changement de direction
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
double ComptageEtat(vector<int> liste, int quest=0, int indiv_focus=0, int etat=1); // Répertorier les états du voisinage d'un individu focus grâce à son etat, sa liste de voisinage (moduler par le cut-off)
double ChooseNewDir(int individu, vector<int> liste);
int Determination_etat(int individu, vector<int> liste);
int Find_head(); // fonction qui trouve l'individu de tête
void Move(double temps, int individu); // Faire bouger les individus actifs
double Dist(int i, int j); // Calculer la distance entre i et j
vector<int> defSousGroupes(int i); // Definir le voisinage de i
void Change_orientation_spont(double temps, double h);
double mod(double val);


//variables
int perc=2, repetition, rep_max=1, rate=100;
double vitesse=0.1, orientation, seuil = 0.01 ;
int N,L, head;
int borne_sup, borne_inf;
double Gam,T,particiones,alpha,beta,Gamma,omega,nss, nos;
double Positions[10000][4], Positions_temp[10000][4], centreGroupei[2],centreGroupe[2];
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
        N           = atoi	(argv[1]);
        T           = strtod(argv[2],&pEnd); //100
        L		= strtod(argv[3],&pEnd); //50
        alpha		= strtod(argv[4],&pEnd); //0.01;
        beta		= strtod(argv[5],&pEnd);//1;
        Gamma		= strtod(argv[6],&pEnd);//0;
        omega		=strtod(argv[7],&pEnd); //0.001;
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

    char filename_Pol[1000];
    strcpy(filename_Pol,"data_Polarization_IGM");
    strcat(filename_Pol,filename);
    ofstream toPol(filename_Pol);
    toPol.precision(5);
    toPol.setf(ios::scientific,ios::floatfield);

    char filename_Coh[1000];
    strcpy(filename_Coh,"data_Cohesion_IGM");
    strcat(filename_Coh,filename);
    ofstream toCoh(filename_Coh);
    toCoh.precision(5);
    toCoh.setf(ios::scientific,ios::floatfield);

    for(repetition=0; repetition<rep_max; repetition++)
    {
        cout << "rep :" << " " << repetition+1 << " /" << rep_max << endl;

        //Las variables de los pasos de tiempo.
        double t = 0.0, h = 0.1;//0.1/pow(Rang,betaCA);

        if(repetition==0) 	// on affiche le h que lors de la première répétition
        {
            cout << "h = " << h << endl;
        }

        //u sirve para calcular un número aleatorio para calcular el paso de tiempo en Gillespie.
        double u = 0.0;

        //N_double sirve para poder dividir cuando se calcula el MSD o la posición del centro de masa.
        double N_double;
        N_double = N;

        Initialisation();
        int nSameState = int(ComptageEtat(defSousGroupes(i))), nDiffState = N-nSameState; // initialisation du compteur d'état

	            centreGroupe[0] =0;
            centreGroupe[1] =0;

        while (t < T)
        {
            double polarizationOrder[2],  polarizationOrderG=0 , cohesionOrder=0;
            polarizationOrder[0]=0;
            polarizationOrder[1]=0;
            centreGroupei[0] =0;
            centreGroupei[1] =0;

            int premier_passage=1;
            borne_inf=0;
            borne_sup=N;

            for(i=borne_inf; i<borne_sup; i++)
            {
                u = ranMT();				//on tire aléatoirement
                nss=ComptageEtat(defSousGroupes(i),1,i);
                nos=ComptageEtat(defSousGroupes(i),1,i,0)-nss;
                //cout << "nss = " << nss << " - nos = " << nos << endl;

                // Update du suivi des configurations
                if(nos==0)  				//tout le monde est dans le même etat que moi
                {
                    if(Positions[i][3]!=1) 			//si mon état précdent n'était pas le même
                    {
                        nDiffState--;
                        nSameState++;
                        Positions[i][3]=1;
                    }
                }
                else 					// Au moins un autre a un état différent
                {
                    if(Positions[i][3]!=0) 		//si mon état précédent n'était pas le même
                    {
                        nDiffState++;
                        nSameState--;
                        Positions[i][3]=0;
                    }
                }


                //Choix de changer ou non de direction
                Gam = (alpha * pow(nos,beta)+omega)/pow(nss,Gamma);
                //cout << Gam << endl;

                if(u < Gam*h)
                {
                    u = ranMT();
                    if(u< (alpha * pow(nos,beta)/(alpha * pow(nos,beta)+omega)) )
                    {
                        Positions[i][2]=ChooseNewDir(i,defSousGroupes(i));
                    }
                    else
                    {
                        Positions[i][2]=ranMT()*2*M_PI;
                    }
                }

                Move(h,i);

            }
            toNM << t << "\t \t" << nSameState << "\t \t" << nDiffState << endl;
            for(i=0; i<N; i++)
            {
                toMovers << t << " " << Positions[i][0] << " " << Positions[i][1] << " " << Positions[i][2] << " " << Positions[i][3] << endl;
                polarizationOrder[0] += cos(Positions[i][2]);	// somme des cos de toutes les orientations des individus
                polarizationOrder[1] += sin(Positions[i][2]);

                if (min(abs(Positions[i][0]-centreGroupe[0]) , min(abs(Positions[i][0]-(centreGroupe[0]-L)), abs(Positions[i][0]-(centreGroupe[0]+L)))) == abs(Positions[i][0]-centreGroupe[0]))
                {
                    centreGroupei[0] += Positions[i][0];
                }
                else if (min(abs(Positions[i][0]-centreGroupe[0]) , min(abs(Positions[i][0]-(centreGroupe[0]-L)), abs(Positions[i][0]-(centreGroupe[0]+L)))) == abs(Positions[i][0]-(centreGroupe[0]-L)))
            	{
                	centreGroupei[0] += Positions[i][0]+L;
                }
                else if (min(abs(Positions[i][0]-centreGroupe[0]) , min(abs(Positions[i][0]-(centreGroupe[0]-L)), abs(Positions[i][0]-(centreGroupe[0]+L)))) == abs(Positions[i][0]-(centreGroupe[0]+L)))
            	{
                	centreGroupei[0] += Positions[i][0]-L;
                }

                if (min(abs(Positions[i][1]-centreGroupe[1]) , min(abs(Positions[i][1]-(centreGroupe[1]-L)), abs(Positions[i][1]-(centreGroupe[1]+L)))) == abs(Positions[i][1]-centreGroupe[1]))
            	{
                	centreGroupei[1] += Positions[i][1];
                }
                else if (min(abs(Positions[i][1]-centreGroupe[1]) , min(abs(Positions[i][1]-(centreGroupe[1]-L)), abs(Positions[i][1]-(centreGroupe[1]+L)))) == abs(Positions[i][1]-(centreGroupe[1]-L)))
            	{
                	centreGroupei[1] += Positions[i][1]+L;
                }
                else if (min(abs(Positions[i][1]-centreGroupe[1]) , min(abs(Positions[i][1]-(centreGroupe[1]-L)), abs(Positions[i][1]-(centreGroupe[1]+L)))) == abs(Positions[i][1]-(centreGroupe[1]+L)))
            	{
                	centreGroupei[1] += Positions[i][1]-L;
                }

                //centreGroupei[0] += Positions[i][0];			// calcul de la position du barycentre du groupe
                                    //centreGroupei[1] += Positions[i][1];

            }
            centreGroupe[0] = mod(centreGroupei[0]/N);
            centreGroupe[1] = mod(centreGroupei[1]/N);
	    //cout << centreGroupe[0] << " " << centreGroupe[1] << endl;
            for(i=0; i<N; i++)
            {
                cohesionOrder+= sqrt(pow( min(abs(Positions[i][1]-centreGroupe[1]) , min(abs(Positions[i][1]-(centreGroupe[1]-L)), abs(Positions[i][1]-(centreGroupe[1]+L)))),2)+pow(min(abs(Positions[i][0]-centreGroupe[0]) , min(abs(Positions[i][0]-(centreGroupe[0]-L)) , abs(Positions[i][0]-(centreGroupe[0]+L))) ),2));; // sommation de toutes les interdistances avec le centre de masse
            }
            polarizationOrderG =  sqrt(pow(polarizationOrder[0],2)+pow(polarizationOrder[1],2))/N;
            toPol << t << " " <<  polarizationOrderG << endl;

            cohesionOrder /= N;
            toCoh << t << " " <<  cohesionOrder << endl;

            t+=h;
            toMovers << endl;
        }

    }
}


vector<int> defSousGroupes(int j)
{
    //cout << "lancement defsousgroupes" << endl;
    vector<int> sousGroupe;
    vector<int>::iterator it;
    sousGroupe.clear();
    int i;
    double tirage;
    for(i=0; i<N; i++) //métrique
    {
        if (i!=j)
        {
            //strong cut-off :
            if( Dist(i,j)<=0.5*perc)
            {
                sousGroupe.push_back(i); //printf("SS grp de %d prend %d\n",j,i);
            }
        }

    }
    //cout << endl;

    return sousGroupe;
}

double ComptageEtat(vector<int> liste, int quest, int indiv_focus, int etat)
{
    //cout << "lancement comptage état" << endl;
    int i;
    double nombreEtat=0;
    if (quest==0)                           //comptage d'état sur la population
    {
        for (i=0; i<N; i++)
        {
            if(Determination_etat(i,liste)==1)
            {
                nombreEtat++;
            }
        }
    }
    else if(quest==1)                      //comptage d'état autour d'un individu
    {
        if (etat==1)
        {
            for (i=0; i<liste.size(); i++)
            {
                if (cos(Positions[indiv_focus][2]-Positions[liste[i]][2])>1-seuil)
                {
                    //cout << cos(Positions[indiv_focus][2]-Positions[liste[i]][2]) << endl;
                    //nombreEtat+=exp(-(Dist(indiv_focus,liste[i]))/(perc));
                    nombreEtat++;
                }
            }
        }
        else
        {
            for (i=0; i<liste.size(); i++)
            {
                nombreEtat++;
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

int Determination_etat(int individu, vector<int> liste)
{
    int etat=1,i;
    for (i=0; i<liste.size(); i++)
    {
        if (cos(Positions[individu][2]-Positions[liste[i]][2])>1-seuil)
        {
            etat=0;
            break;
        }
    }
    return etat;
}


void Move(double temps, int individu)
{
    Positions[individu][0]=mod(Positions[individu][0]+vitesse*cos(Positions[individu][2])*temps);
    Positions[individu][1]=mod(Positions[individu][1]+vitesse*sin(Positions[individu][2])*temps);
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
        dist=sqrt(pow( min(abs(Positions[i][1]-Positions[j][1]) , min(abs(Positions[i][1]-(Positions[j][1]-L)), abs(Positions[i][1]-(Positions[j][1]+L)))),2)+pow(min(abs(Positions[i][0]-Positions[j][0]) , min(abs(Positions[i][0]-(Positions[j][0]-L)) , abs(Positions[i][0]-(Positions[j][0]+L))) ),2));
    }
    //cout << "dist = " << dist <<endl;
    return dist;
}

void Initialisation()
{
    //cout << "lancement initialisation" << endl;
    //Es para utilizar los números aleatorios.
    unsigned int uin = time(NULL);
    srand(uin);
    seedMT(uin);
    int i;
    for (i=0; i<N; i++)
    {
        Positions[i][0]=ranMT()*L;
        Positions[i][1]=ranMT()*L;
        Positions[i][2]=ranMT()*2*M_PI;
        Positions[i][3]=1;
    }
}

double ChooseNewDir(int individu, vector<int> list)
{
    double newdir, prob[N], probTot, tirage;
    int i, j;
    vector<int> liste;
    vector<int>::iterator it;
    liste.clear();

    /*int tirage=ranMT()*list.size();
     newdir=Positions[list[int(tirage)]][2];*/

    for(it=list.begin(); it!=list.end(); ++it) //on remplit la liste de comparaison qui est rechargée avec la liste updatée
    {
        liste.push_back(*it);
    }
    j=0;
    /*for(it=liste.begin(); it!=liste.end(); ++it) //on remplit la liste de comparaison qui est rechargée avec la liste updatée
     {
     cout << *it << " ";
     }
     cout << endl;*/
    for(it=liste.begin(); it!=liste.end(); ++it) // on boucle sur la liste de comparaison
    {
        //cout << *it << " " << (cos(Positions[individu][2]-Positions[*it][2]) < 1-seuil) << endl;
        if( cos(Positions[individu][2]-Positions[*it][2]) > 1-seuil) //si l'individu n'est pas dans le même etat que celui qu'on regarde
        {
            liste.erase(liste.begin()+j);  // on élimine l'individu
            --it;
            j-=1;
        }
        j+=1;
    }

    /*for(it=liste.begin(); it!=liste.end(); ++it)
     {
     cout << *it << " ";
     }*/

    tirage = ranMT(); // on tire respectivement des probabilité de tirage de chacun des individus de la liste
    if (liste.size()!=0)
    {
        for (i=0; i<liste.size()-1; i++)
        {
            prob[i]=1/liste.size();
        }
        j=0;
        probTot=0;
        while (tirage > probTot)
        {
            probTot+=prob[j];
            j+=1;
        }

        newdir=Positions[liste[j-1]][2];

    }
    else
    {
        //cout << "Warning : an individual was asked to change direction without having an heterogeneous neighborhood." << endl;
        newdir=Positions[individu][2];
    }

    return newdir;
}

double mod(double val)  //fonction de traitement des conditions périodiques en position
{
    double valeur;
    if(val<0)
    {
        valeur=val+L;
    }
    else if(val>L)
    {
        valeur=val-L;
    }
    else
    {
        valeur=val;
    }
    return valeur;
}




