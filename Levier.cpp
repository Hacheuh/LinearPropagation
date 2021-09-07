// Ce programme créé une grille, chaque noeud représentant un agent. On génère ensuite leur correspondance topologique puis on créée la matrice d'adjacence associée.
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include <list>
#include <algorithm>
#include <utility>
using namespace std;

extern void seedMT2();
extern void seedMT(unsigned long int);
extern double ranMT(void);

FILE *init, *init_gen, *sortie1, *sortie2, *init_topo, *sortie3;

vector<vector<double> > pos, reseau, topo, mat_dist;    // pos matrice de position, reseau atrice d'adjacence, topo matrice de correspondance topologique
vector<vector<double> >tirage;
vector<double> reseau_l, topo_l;                        // vecteur lignes pour les matrices ci-dessus
int N,percolation, choice, largeur=10;                     // inputs
string Quest="Creation", quest="N";                                     // Variable d'aiguillage : Creation ou lecture de la topologie, Variable d'aiguillage Creation ou lecture de positions.
double vitesse=10.0, delta_t=0.1/vitesse, T=10.0, ran_angle=M_PI/6, D_0, Di_0=1, alpha=0.1;// vitesse des individus mobiles, pas de temps, intervalle temporel, bruit angulaire, constante de distance

double distance(int k, int l, vector<vector<double> > mtx);  // fonction calculant la distance entre deux points de coordonnées xy
void matrice_distance();
double OrderParameter(); 				//Génère le paramètre d'ordre du système

void Conditions_initiales();                        // fonction permettant de choisir entre lecture ou génération
void read();					                    // lire le fichier de position s'il y en a un
void Generate();				                    // générer les positions aléatoirement
void reslist(vector<vector<double> > mtx1);         // fonction générant la matrice la matrice d'adjacence 'reseau' à partir des probabilités d'interaction
void ecriture(double temps);

void Move(double temps, int individu);              // fonction de déplacement des individus mobiles
void Temporality();                                 // fonction lançant le déroulement temporel
double period(double valeur);                       // fonction gérant la périodicité de l'arène

double selecteur_r(int choix, int individu);
double regle_0(int individu);                                     //règle topo-métrique
double regle_1(int individu);                                     //règle Hadrien
double regle_2(int individu);                                     //règle Fernando
double vck(int individu);                                         //règle vicsek

void ChooseCreationOrReading();
void distance_topo();               // Calcul des correspondance topologique : calcul de topo, matrice NxN-1.
void lecture_topo();                // Lecture d'un fichier contenant les correspondance topologique : remplissage de topo

ofstream toRes("resultat_res.txt",ios::out | ios::trunc); // Ouverture des fichiers de sortie
ofstream toPos("resultat_pos.txt",ios::out | ios::trunc);
ofstream toTop("Topologie.txt",ios::out | ios::trunc);
ofstream toOP("OrderParameter.txt",ios::out | ios::trunc);

struct pair_comp // structure définissant la règle de comparaison du sort pour des <pair>
{
    int operator()(const pair<int,double>& t1,const pair<int,double>& t2) const
    {
        return t1.second<t2.second;
    }
};

int main(int argc, const char * argv[])
{
    int i,j,k;

    if (argc > 1)// 2 paramètres : N la taille de pop, le coefficient de percolation : nombre de voisins pris en compte
    {
        char * pEnd;
        N               = strtod(argv[1],&pEnd);
        percolation		= strtod(argv[2],&pEnd);
        D_0             = strtod(argv[3],&pEnd);
        choice		    = strtod(argv[4],&pEnd);

        //Quest           = strtod(argv[5],&pEnd); //Paramètre optionnel d'aiguillage si l'on veut switch souvent entre les deux possibilités
    }

    unsigned int uin = time(NULL);
    srand(uin);
    seedMT(uin);

    Conditions_initiales();     // Génération des Coordonnées spatiales des noeuds de la grille
    Temporality();

    toRes.close();              //Fermeture des fichiers de sortie
    toPos.close();

}

void ChooseCreationOrReading()
{
    if (Quest=="Creation")
    {
        distance_topo();
    }
    else if (Quest=="Reading")
    {
        lecture_topo();
    }
    else
    {
        cout << "You didn't choose any of the two proposed choices. Please choose Reading mode by typing /'Reading/' or Creation mode by typing /'Creation/'." << endl;
    }
}

void Conditions_initiales()
{
    if (quest=="N") //on ne veut pas générer de nouvelles conditions initiales, on lance la lecture du fichier par defaut
    {
        read();
    }
    else if (quest=="Y")
    {
        Generate();//on veut générer de nouvelles conditions initiales, on lance la fonction generate
    }
    else
    {
        printf("You didn't type Y or N\n");
        exit(1); // si aucune des deux solutions n'a été tapé le script s'arrête
    }
}


void Generate()
{
    int i, j;
    vector<double> pos_l;
    init_gen=fopen("initiation_gen.txt","w+");
    for (i=0; i<N; i++)
    {
        for (j=0; j<2; j++)
        {
            pos_l.push_back(ranMT()*10-5);                        // On remplit pour chaque individu son x et son y
        }
        pos_l.push_back(ranMT()*2*M_PI-M_PI);                          // On génère un angle
        pos.push_back(pos_l);                                   // On stocke ce vecteur dans le vecteur pos
        pos_l.clear();                                          // On vide pos_l
        fprintf(init_gen, "%lf %lf %lf\n",pos[i][0],pos[i][1], pos[i][2]);    // On écrit la ligne dans un fichier
    }
    fclose(init_gen);
}

void read()
{
    int i,j;
    double stock;
    vector<double> pos_l;
    init=fopen("initiation_gen.txt","r");
    for(i=0; i<N; i++)
    {
        for(j=0; j<3; j++)
        {
            fscanf(init,"%lf", &stock); //On lit le fichier nombre par nombre en stockant chaque nombre dans une variable
            pos_l.push_back(stock);     //cette variable est ajoutée au vecteur ligne pos_l
        }
        pos.push_back(pos_l);           //Tous les trois nombres, on ajoute le tout à pos
        pos_l.clear();                  //On vide pos_l
    }
    fclose(init);
}


double distance(int k, int l, vector<vector<double> > mtx)
{
    double dist_x, dist_y;
    dist_x=min(sqrt(pow(mtx[k][0]-mtx[l][0],2)),sqrt(pow(mtx[k][0]-mtx[l][0]+2*largeur,2)));
    dist_y=min(sqrt(pow(mtx[k][1]-mtx[l][1],2)),sqrt(pow(mtx[k][1]-mtx[l][1]+2*largeur,2)));

    return sqrt(pow(dist_x,2)+pow(dist_y,2)); // calcul de distance métrique entre deux points
}

void reslist(vector<vector<double> > mtx1)
{
    int i,j,k,tir2,relancer;
    int temp;
    double tir;
    vector<double>::iterator it, it2;
    vector<vector<double> >::iterator it3;
    reseau.clear();

    vector<double> tirage_l, possibilities;
    // Remplissage de reseau pour obtenir la taille NxN
    for(i=0; i<N; i++)
    {
        for(j=0; j<N; j++)
        {
            reseau_l.push_back(0);
        }
        reseau.push_back(reseau_l);
    }

    //Creation des liens dans la matrice d'adjacence 'reseau'
    if(choice==2)
    {
        tirage.clear(); // on refait entièrement le vecteur de voisinage aléatoire dans la condition 2 uniquement
    }

    for(i=0; i<N; i++)
    {
        if(choice==3)                                       //dans le cas où on veut voir le vicsek pur le reseau est constitué différemment
        {
            j=0;
            while(mat_dist[i][topo[i][j]]<D_0)
            {
                temp = topo[i][j];              // On récupère l'ID de l'agent
                reseau[i][temp]=1;              // On crée le lien entre i et temp dans la matrice d'adjacence
                j++;
            }
        }
        else if(choice==2)
        {

            tirage_l.clear();
            possibilities.clear();
            for(j=0; j<N; j++)
            {
                if(j!=i)
                {
                    possibilities.push_back(j);
                }
            }
            random_shuffle ( possibilities.begin(), possibilities.end() );
            for(j=0; j<percolation; j++)
            {
                tirage_l.push_back(possibilities[j]);
                reseau[i][possibilities[j]]=1;
            }
            tirage.push_back(tirage_l);
        }
        else if(choice==4) // dans la regle 4, on tire aléatoirement selon alpha le renouvellement ou non du voisinage
        {
            int check;
            tirage_l.clear();
            possibilities.clear();

            for(j=0; j<percolation; j++)
            {
                check=0;                // condition de verification qu'on ait pas repris le même ou tiré parmis les autres déjà dans le vecteur
                tir=ranMT();            //tirage du renouvellemnt
                if(tir<alpha && percolation<N-1)           // verification de la condition de renouvellement, si percolation =N-1 boucle infinie car personne ne peut être pris pour remplacer
                {
                    while (check<percolation) // tant qu'on a pas trouvé un bon remplaçant
                    {
                        tir2=ranMT()*N;     // on tire un nouvel individu
                        for(k=0; k<percolation; k++)  // on parcours l'ancien vecteur de voisin tiré
                        {
                            if(tir2!=tirage[i][k] && tir2!=i)      //s'il n'y a pas égalité avec le nouveau tiré ...
                            {
                                check+=1;               // check s'incrémente
                            }

                        }
                        if(check==percolation)          //si autant de check qu'on a de voisin possible
                        {
                            tirage[i][j]=tir2;          // alors on update ce voisin avec le nouveau tiré
                            reseau[i][tir2]=1;
                        }
                        else                            // sinon check se remet à 0 et on retire
                        {
                            check=0;
                        }
                    }
                }
                else{
                    reseau[i][tirage[i][j]]=1;
                }
            }
        }
        else
        {
            for(j=0; j<floor(percolation); j++) // On tronque la liste de aux 'percolation' premiers agents
            {
                temp = topo[i][j];              // On récupère l'ID de l'agent
                reseau[i][temp]=1;              // On crée le lien entre i et temp dans la matrice d'adjacence
            }
        }
    }
    /*cout << "tirage" << endl;
    for(it3=tirage.begin(); it3!=tirage.end(); ++it3)
    {
        for(it=it3->begin(); it!=it3->end(); ++it)
        {
            cout << *it << " ";
        }
        cout << endl;
    }
    cout << endl;*/
}

void ecriture(double temps)
{
    int i,j;
    toRes << temps << endl;
    toPos << temps << endl;
    for (i=0; i<N; i++)
    {
        for (j=0; j<N; j++)
        {
            //fprintf(sortie1, "%lf %lf %lf \n", reseau[i][j], distance(i,j,pos), proba_interact(i,j));
            toRes << reseau[i][j] << " ";
        }
        toRes << endl;
        toPos << pos[i][0] << " " << pos[i][1] << " " << pos[i][2] << endl;
    }

    toOP << OrderParameter() << endl;
    toRes << endl;
    toPos << endl;

}

void distance_topo()
{
    pair_comp comp;     // foncteur règle de comparation
    int i,j,k,l;
    vector<double> topolo_l, sorted;
    vector<pair<int, double> > unsorted;
    vector<pair<int, double> >::iterator it;
    topo.clear();

    for(i=0; i<N; i++)                                              // On parcourt chaque agent...
    {
        unsorted.clear();                                           // On (ré)initialise le vecteur de pairs
        for(j=0; j<N; j++)                                          // ... et on regarde chaque individus ...
        {
            if (i!=j)                                               // ... qui n'est pas soi-même ...
            {
                unsorted.push_back(make_pair(j,mat_dist[i][j]));    // ... on crée la paire individu/distance associée, et on l'ajoute au vecteur de paire de i ...
            }
        }
        sort(unsorted.begin(),unsorted.end(),comp);                 // ... on trie le vecteur de paire suivant sa deuxième colonne (règle comp) ...
        sorted.clear();                                             // On (ré)initialise le vecteur de liste topologique triée
        for(it=unsorted.begin(); it!=unsorted.end(); ++it)          // ... on rempli le vecteur topologique de i avec la deuxième colonne de son vecteur de paire ...
        {
            sorted.push_back(it->first);
        }
        topo.push_back(sorted);                                     // ... on push ce vecteur ligne dans la matrice topologique.
    }

    // Ecriture du vecteur topologique dans son fichier.
    for(i=0; i<N; i++)
    {
        for(j=0; j<N-1; j++)
        {
            toTop << topo[i][j] << " ";
            //cout << topo[i][j] << " ";
        }
        toTop << endl;
        //cout << endl;
    }
    toTop << endl;
}

void lecture_topo()
{
    int i,j;
    ifstream fromTop("Topologie.txt",ios::in);

    // remplissage de la matrice topologique pour obtenir la taille NxN-1
    for(i=0; i<N; i++)
    {
        topo_l.clear();
        for(j=0; j<N-1; j++)
        {
            topo_l.push_back(0);
        }
        topo.push_back(topo_l);
    }

    // Si le fichier ouvre, remplissage de la matrice topologique
    if(fromTop)
    {
        for (i=0; i<N; i++)
        {
            for (j=0; j<N-1; j++)
            {
                fromTop >> topo[i][j];
            }
        }
    }
}


void matrice_distance()
{
    int i,j;
    vector<double> mat_dist_l;
    mat_dist.clear();
    // remplissage de la matrice de distance pour obtenir la taille NxN
    for(i=0; i<N; i++)
    {
        for(j=0; j<N; j++)
        {
            mat_dist_l.push_back(0);
        }
        mat_dist.push_back(mat_dist_l);
    }

    // On calcul la distance pour chaque (i,j) // Optimisation possible car symétrie
    for(i=0; i<N; i++)
    {
        for(j=i; j<N; j++)
        {
            if (i==j)
            {
                mat_dist[i][j]=0;
            }
            else
            {
                mat_dist[i][j]=distance(i,j,pos);
                mat_dist[j][i]=distance(i,j,pos);
                //cout << mat_dist[i][j] << " ";
            }
        }
        //cout << endl;
    }
}


void Move(double temps, int individu)
{
    pos[individu][0]=period(pos[individu][0]+vitesse*cos(pos[individu][2])*temps);
    pos[individu][1]=period(pos[individu][1]+vitesse*sin(pos[individu][2])*temps);
    pos[individu][2]+=(selecteur_r(choice,individu)*temps);//+(sqrt(2*Di_0)*(ranMT()*2*M_PI)));

}

void Temporality()
{
    int i,c_temp;
    double t=0.0;
    while (t<T)
    {
        matrice_distance();         // Creation de la matrice contenant toutes les distances inter-individuelles
        if(choice!=2 || choice!=4)
        {
            ChooseCreationOrReading();  // Aiguillage en fonction de la variable Quest
        }
        if(t==0 && choice==4)
        {
            c_temp=choice;
            choice=2;
        }
        reslist(reseau);            // Construction du reseau

        for(i=0; i<N; i++)
        {
            Move(delta_t,i);
        }
        if(t==0 && choice!=c_temp)
        {
            choice=c_temp;
        }

        ecriture(t);                 // Ecriture des fichiers positions et reseau.
        t+=delta_t;
    }
}

double period(double valeur)
{
    double periodic;
    if (valeur>=largeur)
    {
        periodic=valeur-2*largeur;
    }
    else if(valeur<=-largeur)
    {
        periodic=valeur+2*largeur;
    }
    else
    {
        periodic=valeur;
    }
    return periodic;
}

double selecteur_r(int choix, int individu)
{
    double selection;
    switch (choix)
    {
    case 0 :
        selection=regle_0(individu);
        break;
    case 1 :
        selection=regle_1(individu);
        break;
    case 2 :
        selection=regle_2(individu);
        break;
    case 3 :
        selection=vck(individu);
        break;
    case 4 :
        selection=regle_2(individu); // même règle mais voisinage choisi différement
        break;
    default :
        cout << "Input choix invalide" << endl;
        break;
    }
    return selection;
}

double regle_0(int individu)
{
    int i;
    double angle_moyen=0;
    for(i=0; i<floor(percolation); i++)
    {
        angle_moyen+=sin(pos[topo[individu][i]][2]-pos[individu][2]);
    }
    //angle_moyen/=somm;
    return angle_moyen;
}

double regle_1(int individu)
{

}

double regle_2(int individu)
{
    int i,j;
    double somm=0, angle_moyen=0;
    for(j=0; j<floor(percolation); j++)
    {
        angle_moyen+=(sin(pos[tirage[individu][j]][2]-pos[individu][2])*exp(-mat_dist[individu][tirage[individu][j]]/D_0));
        somm+=exp(-mat_dist[individu][tirage[individu][j]]/D_0);
    }
    angle_moyen/=somm;
    return angle_moyen;
}

double vck(int individu)
{
    int i=0;
    double angle_moyen=0;
    while(mat_dist[individu][topo[individu][i]]<D_0)
    {
        angle_moyen+=sin(pos[topo[individu][i]][2]-pos[individu][2]);
        i++;
    }
    if (i==0)
    {
        angle_moyen=0;
    }
    else
    {
        angle_moyen/=i;
    }
    return angle_moyen;
}

double OrderParameter()
{
    int i;
    double OrP,OrPx=0,OrPy=0;
    for(i=0; i<N; i++)
    {
        OrPx+=cos(pos[i][2]);
        OrPy+=sin(pos[i][2]);
    }
    OrP=sqrt(pow(OrPx/N,2)+pow(OrPy/N,2));
    return OrP;
}
