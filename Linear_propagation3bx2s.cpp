// On regarde avec ce programme la propagation d'un état le long d'une ligne d'agents. La version 2 introduit un état mobile et donc un tracking dans l'espace des dernières particules en mouvement. La version 3, modifie fondamentalement la façon de prévoir les évènements, on intégre de manière brut le système temporellement.



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

int N,largeur,a;

const int nombreEvenement=2000;
int evenementControl=1, tester=0; // à combien on commence,
double Rang =1.53;
int perc=2, repetition, rep_max=1;
double vitesse=-0.5, perc_moyen=0.0; // nombre moyen d'individu perçu
vector<double> perc_t; // vecteur nous permettant d'update le h

int ActivationMovers=0; //1 pour activer l'écriture du fichier movers

vector<int> defSousGroupes(int i); // Definir le voisinage de i
vector<int> Etats, Etats_temp;
vector<double> gamma_indiv(double alphaAB, double betaAB, double gammaAB, double omegaAB, double alphaBC, double betaBC, double gammaBC, double omegaBC, double alphaCA, double betaCA, double gammaCA, double omegaCA,  double N, vector<int> etat); // Definir tous les taux de transitions par unité de temps pour les N individus avec les paramètres du modèle de transition (alpha beta gamma..) et l'etat du système
vector<int> Liste_etat(), Movers;
double ComptageEtat(int etat, vector<int> liste, int quest=0, int indiv_focus=0); // Répertorier les états du voisinage d'un individu focus grâce à son etat, sa liste de voisinage (moduler par le cut-off)
double ProbAB(double alpha, double beta, double gamma, double omega, double N, int i); // Calculer la proba de transition AB
double ProbBC(double alpha, double beta, double gamma, double omega, double N, int i); // Calculer la proba de transition BC
double ProbCA(double alpha, double beta, double gamma, double omega, double N, int i); // Calculer la proba de transition CA
double Positions[10000], Positions_temp[10000];
void Move(double temps, int individu); // Faire bouger les individus actifs
vector<vector<double> > MoversPos;
vector<double> MoversPos_temp;
double Dist(int i, int j); // Calculer la distance entre i et j
double findMax(vector<double> vec); // Trouve le maximum des valeurs d'un vecteur



//Funciones externas a main() para calcular los rates.
double ProbAB(double alpha, double beta, double gamma, double omega, double N, int i)
{
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
{
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
{
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


//Se leen los parámetros en la línea de comando.
int main(int argc, const char * argv[])
{
    //Es para utilizar los números aleatorios.
    unsigned int uin = time(NULL);
    srand(uin);
    seedMT(uin);

    int i,j,k,l;
    int front_mover, last_front_mover, borne_sup, borne_inf, first_front_static; // Definition de point de repère spatiaux (le premier statique, le dernier en mouvement, les bornes d'optimisation)
    //int N;
    double T,particiones,omegaAB,omegaBC,omegaCA,alphaAB,alphaBC,alphaCA,betaAB,betaBC,betaCA,gammaAB,gammaBC,gammaCA;
    double delta; // variable d'écart entre premier statique et dernier en mouvement

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
    strcpy(filename_nm,"data_nM_lipro");
    strcat(filename_nm,filename);
    ofstream toNM(filename_nm);
    toNM.precision(5);
    toNM.setf(ios::scientific,ios::floatfield);

    char filename_FrontVague[1000];
    strcpy(filename_FrontVague,"data_FrontVague_lipro");
    strcat(filename_FrontVague,filename);
    ofstream toFrontVague(filename_FrontVague);
    toFrontVague.precision(5);
    toFrontVague.setf(ios::scientific,ios::floatfield);

    char filename_FrontVagueL[1000];
    strcpy(filename_FrontVagueL,"data_FrontVagueLeft_lipro");
    strcat(filename_FrontVagueL,filename);
    ofstream toFrontVagueLeft(filename_FrontVagueL);
    toFrontVagueLeft.precision(5);
    toFrontVagueLeft.setf(ios::scientific,ios::floatfield);

    char filename_Movers[1000];
    strcpy(filename_Movers,"data_Movers_lipro");
    strcat(filename_Movers,filename);
    ofstream toMovers(filename_Movers);
    toMovers.precision(5);
    toMovers.setf(ios::scientific,ios::floatfield);

    char filename_Delta[1000];
    strcpy(filename_Delta,"data_Delta_lipro");
    strcat(filename_Delta,filename);
    ofstream toDelta(filename_Delta);
    toDelta.precision(5);
    toDelta.setf(ios::scientific,ios::floatfield);

    char filename_Finalid[1000];
    strcpy(filename_Finalid,"data_Finalid_lipro");
    strcat(filename_Finalid,filename);
    ofstream toFinal(filename_Finalid);
    toFinal.precision(5);
    toFinal.setf(ios::scientific,ios::floatfield);

    char filename_Mass[1000];
    strcpy(filename_Mass,"data_Mass_lipro");
    strcat(filename_Mass,filename);
    ofstream toMass(filename_Mass);
    toMass.precision(5);
    toMass.setf(ios::scientific,ios::floatfield);


    for(repetition=0; repetition<rep_max; repetition++)
    {
        cout << "rep :" << " " << repetition+1 << " /" << rep_max << endl;

        //Las variables de los pasos de tiempo.
        double t = 0.0, h = 0.1/pow(Rang,betaCA);

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

        //Variables para calcular los números de ocupación.
        Etats.clear();
        Etats=Liste_etat(); // initialisation des etats
        /*for(i=0; i<N; i++)
        {
            printf("%d ", Etats[i]);
        }
        printf("\n");*/

        for (i=0; i<N; i++) // initialisation des positions
        {
            Positions[i]=i;
        }

        Movers.clear();
        MoversPos.clear();
        Movers.push_back(0); // La liste des individus en déplacement contient au moins le premier i.e. 0
        last_front_mover=N/2; // Le dernier en mouvement de la ligne est 0 au départ
        front_mover=N/2; // même raison
        first_front_static=N/2+1; // Le premier statique est forcément 1

        int last_front_moverL=N/2; // Le dernier en mouvement de la ligne est 0 au départ côté gauche
        int front_moverL=N/2; // même raison
        int first_front_staticL=N/2-1; // Le premier statique à gauche est forcément -1

        evenementControl=1; // On initialise le nombre d'evenement à 1, la particule 0 a virtuellement déjà eut sa transition
        vector<int>::iterator it;
        vector<double>::iterator itd;

        perc_moyen=0.0;

        int nA = int(ComptageEtat(0, Etats)), nB = int(ComptageEtat(1, Etats)), nC = int(ComptageEtat(2, Etats)); // initialisation du compteur d'état
        //printf("%d %d %d\n", nA, nB, nC);

        while (t < T)
        {

            int premier_passage=1; //premier passage dans le sens gauche droite
	    int premier_passageL=1; // premier passage dans le sens droite gauche
	    perc_t.clear();
            vector<double> gammas = gamma_indiv(alphaAB,betaAB,gammaAB,omegaAB,alphaBC,betaBC,gammaBC,omegaBC,alphaCA,betaCA,gammaCA,omegaCA,N_double,Etats); // on génére les taux de transitions de toutes les particules à cet instant. Les distances sont calculer dans cet appel.


            /*if(last_front_mover+5*perc>=N) 		// Optimisation de boucle pour pas prendre en compte tout ceux qui ne peuvent pas bouger.
            {
                borne_sup=N;
            }
            else
            {
                borne_sup=last_front_mover+5*perc;
            }

            if(last_front_mover-5*perc<=0) 		// Optimisation de boucle pour pas prendre en compte tout ceux qui ont bougé les premiers mais qui n'apporte rien à la suite.
            {
                borne_inf=0;
            }
            else
            {
                borne_inf=last_front_mover-3*perc;
            }*/

            borne_inf=0;
            borne_sup=N;

            for(i=borne_inf; i<borne_sup; i++)
            {
                if(Etats[i]==0)
                {
                    Move(h,i);			//si c'est un état mobile, on intégre le déplacement
                }
                else if(Etats[i]==2)
                {
                    u = ranMT();			//on tire aléatoirement
                    //if (i<20){cout << i << "\t \t" << u << " " << gammas[i]*h << "\t \t" << gammas[i] << "\t \t" << endl;}
                    if(u < gammas[i]*h) 		//si le tirage est OK
                    {
                        //cout << (gammas[i]*h) << endl;
                        Etats[i]=0;		//l'individu passe dans l'état mobile
                        nA++;
                        nC--;
                        Movers.push_back(i);	//on l'ajoute à la liste des mobiles
                        front_mover=i;		//le dernier i à bouger(puisqu'on parcours dans le sens croissant) est le front_mover de cette vague d'événement
			if(premier_passageL=1){
				front_moverL=i; 	// on attrape le premier qui declenche la transition dans l'ordre de parcours 0 -> N-1
				premier_passageL=0;
			}
                        evenementControl++;	//on ajoute autant d'événement qu'il y a de transitions d'état
                        perc_moyen+=ComptageEtat(2, defSousGroupes(i),1,i);
			perc_t.push_back(ComptageEtat(2, defSousGroupes(i),1,i));

                    }
                    else if(premier_passage==1)//on attrappe le premier qui n'a pas declenché la transition
                    {
                        //cout << i << "\t \t" << last_front_mover << "\t \t" << (gammas[i]*h) << endl;
                        first_front_static=i;
                        premier_passage=0;
                    }
                }
                else
                {
                    cout << "Anomalie d'état détectée" << endl;
                    return 0;
                }
            }
            //cout << first_front_static << " " << (gammas[first_front_static]*h) << " " << Dist(first_front_static, first_front_static-1) << " " << exp(-(Dist(first_front_static, first_front_static-1))/(perc)) << endl;
            if(last_front_mover < front_mover)  	//on regarde si le front_mover de la vague dépasse les front_movers des autres vagues
            {
                last_front_mover=front_mover;
            }
            if(last_front_moverL > front_moverL)  	//on regarde si le front_mover de la vague dépasse les front_movers des autres vagues
            {
                last_front_moverL=front_moverL;
            }
            if (ActivationMovers==1)
            {
                MoversPos_temp.clear();
                for(it=Movers.begin(); it!=Movers.end(); ++it)
                {
                    MoversPos_temp.push_back(Positions[*it]);
                }
                MoversPos.push_back(MoversPos_temp);
            }
            toNM << t << "\t \t" << nA << "\t \t" << nB << "\t \t" << nC << endl;

            toFrontVague << t << "\t \t";
            toFrontVague << last_front_mover << endl;

            toFrontVagueLeft << t << "\t \t";
            toFrontVagueLeft << last_front_moverL << endl;

            delta=-Positions[last_front_mover]+Positions[first_front_static];
            toDelta << t << "\t \t";
            toDelta << delta << endl;

            if (evenementControl>=nombreEvenement) // pour terminer la simulations si on a fini de propager à toute la ligne
            {
                printf("Dépassement de l'événement contrôle\n");
                break;
            }

            int compteur_perc=0;
            for(i=0; i<N; i++)
            {
                if (gammas[i]>=0.001)
                {
                    compteur_perc=1;
                }
            }
            if(compteur_perc<1)//si au moins un individu a une proba supérieure à 1/100 par pas de temps de changer d'état alors on continue
            {
                printf("Dépassement de perception : %d\n", first_front_static);
                break;
            }

	    if (findMax(perc_t)!=0){
	    	h=0.1/pow(max(findMax(perc_t),1.53),betaCA);
	    }
            t+=h;

        }
        perc_moyen/=evenementControl;
        printf("Le nombre moyen d'individu perçu est de : %lf\n", perc_moyen);
        printf("La durée de l'xp est de : %lf\n",t);
        if(rep_max!=1) // saut de ligne entre chaque répétition
        {
            toDelta << endl;
            toNM << endl;
            toFrontVague << endl;
        }
        toFinal << first_front_static << endl;
	toMass << nA << endl;

        if (ActivationMovers==1)
        {
            for (i=0; i<MoversPos.size(); i++)
            {
                for(itd=MoversPos[i].begin(); itd!=MoversPos[i].end(); ++itd)
                {
                    toMovers << *itd << "\t \t";
                }
                toMovers << endl;
            }
            toMovers << endl;
        }

    }

    toNM.close();

    cout << "Fin du programme." << endl;
    return 0;
}

vector<int> defSousGroupes(int j)
{
    vector<int> sousGroupe;
    sousGroupe.clear();
    int i;
    double tirage;
    //periodicité horizontale topologique
    /*for(i=1; i<1+1; i++) // Attention ce rang doit être un int
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

vector<double> gamma_indiv(double alphaAB, double betaAB, double gammaAB, double omegaAB, double alphaBC, double betaBC, double gammaBC, double omegaBC, double alphaCA, double betaCA, double gammaCA, double omegaCA,  double N, vector<int> etat)
{
    int i;
    vector<double> gammas;
    gammas.clear();
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
        if(i==N/2)
        {
            Liste.push_back(0); //Le premier est actif
        }
        else
        {
            Liste.push_back(2); //Les autres sont inactifs
        }
    }
    return Liste;
}

double ComptageEtat(int etat, vector<int> liste, int quest, int indiv_focus)
{
    int i;
    double nombreEtat=0;
    double tirage;
    if (quest==0)                           //calcul du nombre d'état de manière default, pour l'initialisation des nM
    {
        for (i=0; i<liste.size(); i++)
        {

            if (liste[i]==etat)
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
            if (Etats[liste[i]]==etat)
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
    vector<int>::iterator it;

    Positions[individu]+=(vitesse*temps);
    //cout << Positions << endl;
    //cout << Positions[0] << endl;
}

double Dist(int i, int j)
{
    double dist;
    dist=sqrt(pow(Positions[i]-Positions[j],2));
    return dist;
}


double findMax(vector<double> vec)
{
  double maximum =0;
  vector<double>::iterator it;
  for(it=vec.begin();it!=vec.end();++it)
    {
      if (*it>maximum)
    	maximum = *it ;
    }

  return maximum;
}





